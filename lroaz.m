function [pWorst,x0,lambda0,mu0,indivScens,zupper] = lroaz_version13(gammaprime, nfactor, periods1)

clear get_stage_vectors

if nargin < 3
    periods1 = 10;
    if nargin < 2
        nfactor = 5;
        if nargin < 1
            gammaprime = 0.19;
        end
    end
end

% lroaz is an attempt to make an LRO for the model of Tucson that I have.
% It is based on lroslp_version12.
% I have added some input and output arguments to this version

% optimizer = 'fmincon';
optimizer = 'linprog';
% optimizer = 'cvx';

ConnectionsFile = 'all_scenarios/5/Connections.xlsx';
cellInputFile = { ...
    'all_scenarios/5/Inputs.xlsx', ...
    'all_scenarios/6/Inputs.xlsx', ...
    'all_scenarios/7/Inputs.xlsx', ...
    'all_scenarios/8/Inputs.xlsx', ...
    };

% Problem Parameters
numscen = nfactor*ones(size(cellInputFile));
N = sum(numscen);
% The constant inside the log is a number less than 1
% gammaprime = 0.95;
Nbar = N*(log(N)-1) - log(gammaprime);
Period = 41;
% periods1 = 10;
% c = 1;

[c,A,Rhs,l,u] = get_stage_vectors(1,1, ...
    ConnectionsFile,cellInputFile,Period,periods1);
[x0 initCost] = linprog( c, [], [], A, Rhs, l, u );
A = [A zeros(size(A,1),3)];

% Initialize everything
% x0 = -1;
lambda0 = 1;
zlower = -Inf;
zupper = Inf;
objA = [];
objRhs = [];
feasA = [];
feasRhs = [];
feasSlope = [];
feasInt = [];

switch optimizer
    case 'fmincon'
        options = optimset('MaxIter',1000, ...
            'Algorithm','interior-point', ...
            'GradObj','on');
    case 'linprog'
        options = optimset('MaxIter',85);
    case 'cvx'
        
    otherwise
        error(['Optimizer ' optimizer ' not supported'])
end

% Uniform lower bounds on scenario costs
scenLowBnd = 0;
tic

% Solve all the subproblems & find a good starting mu
[indivScens slope intercept] = solve_scens(x0,numscen);

% Rescaling code:
scale = rescale_problem(initCost,indivScens,numscen,Nbar);
get_stage_vectors('scale',scale);
c = get_stage_vectors(1);
indivScens = scale*indivScens;

mu0 = find_mu(lambda0, numscen,indivScens);

% Initialize variables for exit conditions
pWorst = lambda0*numscen./(mu0-indivScens');
tolerance = 1e-5;
notPrimalSolnFound = true;

% Get the initial objective cut
[theta0 slope intercept] = get_cut(x0,lambda0,mu0,numscen,N,Nbar,indivScens,slope,intercept);

% Variables for trust region
trustRegion = max(abs([x0;lambda0;mu0;theta0]));
rhoBound = 1/4;
scaleDown = 1/4;
scaleUp = 3;
eta = 1/5;
trustRatio = 0.99;

while notPrimalSolnFound || abs(1-sum(pWorst)) > 1e-3
%     Update the matrices of objective and feasibility cuts
    objA = [objA; slope, -1];
    objRhs = [objRhs; -intercept];
    feasA = [feasA; feasSlope];
    feasRhs = [feasRhs; -feasInt];
    
    

%     Format for decision variables: [x lambda mu theta]
    initGuess = [x0; lambda0; mu0; theta0];
    xOld = x0;
    lambdaOld = lambda0;
    muOld = mu0;
    thetaOld = theta0;
    
%     Define bounds
    lowerTrust = initGuess - trustRegion;
    lowerTrust(end) = -Inf;
    upperTrust = initGuess + trustRegion;
    upperTrust(end) = Inf;
    lowerBound = [l;0;scenLowBnd;-Inf];
%     upperBound = [u;max(1,10*lambda0);10*mu0;Inf];
    upperBound = [u;Inf;Inf;Inf];
    lowerBound = max(lowerBound, lowerTrust);
    upperBound = min(upperBound, upperTrust);
    
    disp([ num2str(size(objA,1)) ' objective cuts, ' ...
        num2str(size(feasA,1)) ' feasibility cuts'])
    toc
%     Solve the master problem and get the solution
    switch optimizer
        case 'fmincon'
            [x,~,exitflag] = fmincon( @(x) opt_obj(x,c,N,Nbar), initGuess, ...
                [objA; feasA], [objRhs; feasRhs], A, Rhs, ...
                lowerBound, upperBound, [], options );
        case 'linprog'
            linobj = linear_obj(c,Nbar);
            [x,~,exitflag] = linprog( linobj, ...
                [objA; feasA], [objRhs; feasRhs], A, Rhs, ...
                lowerBound, upperBound, ...
                initGuess, options);
        case 'cvx'
            linobj = linear_obj(c,Nbar);
            cvx_begin quiet
                variable x(length(initGuess))
                cvx_precision best
                minimize linobj*x
                subject to
                    [objA; feasA]*x <= [objRhs; feasRhs];
                    A*x == Rhs;
                    lbIndex = find(lowerBound > -Inf);
                    ubIndex = find(upperBound < Inf);
                    lowerBound(lbIndex) <= x(lbIndex);
                    x(ubIndex) <= upperBound(ubIndex);
%                     lowerBound <= x <= upperBound;
            cvx_end
            exitflag = 1*strcmp(cvx_status,'Solved') + ...
                -2 * strcmp(cvx_status,'Infeasible') + ...
                -3 * strcmp(cvx_status,'Unbounded' ) + ...
                 6 * strcmp(cvx_status,'Inaccurate/Solved');
            exitflag = exitflag - 40*(exitflag==0);
        otherwise
            error(['Optimizer ' optimizer ' not supported'])
    end
    x0 = x(1:end-3);
    lambda0 = x(end-2);
    mu0 = x(end-1);
    thetaMaster = x(end);
    disp(['Scenario Observations: ' num2str(numscen)])
    disp(['Exit flag is ' num2str(exitflag) '.'])
    
%     Solve Subproblems
    [indivScens slope intercept] = solve_scens(x0,numscen);
    
%     If mu infeasible, generate feasibility cut and prevent zlower from
%     updating
    [hMax hIndex] = max(indivScens);
    muFeasible = mu0 > hMax;
    if ~muFeasible
        %         If mu infeasible, generate feasibility cuts and find feasible mu
        [feasSlope feasInt] = get_feas_cut(slope,intercept,hIndex);
        mu0 = find_mu(lambda0, numscen,indivScens);
        disp('Generating feasibility cut.')
        %         plot_feas_step(feasSlope,-feasInt,hIndex);
    else
        feasSlope = [];
        feasInt = [];
    end

%     Get true second stage cost and a new cut
    [theta0 slope intercept] = get_cut(x0,lambda0,mu0,numscen,N,Nbar,indivScens,slope,intercept);
    
%     Calculate whether inside the interior of the trust region
%     trustRegionInterior = (lambda0 < 0.9*upperBound(end-2)) & (mu0 < 0.9*(upperBound(end-1)));
    upperT = ((1+trustRatio)*upperTrust(1:end-1) + (1-trustRatio)*lowerTrust(1:end-1)) / 2;
    lowerT = ((1-trustRatio)*upperTrust(1:end-1) + (1+trustRatio)*lowerTrust(1:end-1)) / 2;
    trustRegionInterior = ~(sum(lowerT > x(1:end-1)) || sum(x(1:end-1) > upperT));
    
%     Modify trust region size
    rho = (thetaOld - theta0) / (thetaOld - thetaMaster);
%     assert(rho <= 1, ['rho = ' num2str(rho) ' > 1, thetaMaster = '...
%         num2str(thetaMaster) ', theta0 = ' num2str(theta0)])
    if rho < rhoBound
        trustRegion = scaleDown * trustRegion;
        disp(['Trust region scaled down to ' num2str(trustRegion)])
    elseif rho > 1-rhoBound && ~trustRegionInterior
        trustRegion = scaleUp * trustRegion;
        disp(['Trust region scaled up to ' num2str(trustRegion)])
    end
    if rho > eta
        newSolution = true;
        disp('New solution found')
    else
        newSolution = false;
        x0 = xOld;
        lambda0 = lambdaOld;
        mu0 = muOld;
        theta0 = thetaOld;
    end
    
%     Update zlower if successful optimization in the interior of the trust
%     region
    switch exitflag
        case 1
            if muFeasible && trustRegionInterior && newSolution
                zlower = get_first_stage_obj(x0,lambda0,mu0,c,N,Nbar)+thetaMaster;
                disp(['Feasible solution, updating zlower = ' num2str(zlower)])
            end
        case 0
            switch optimizer
                case {'linprog','fmincon'}
                    options = optimset(options,'MaxIter',2*options.MaxIter);
                    disp(['Number of iterations increased to ' num2str(options.MaxIter) '.'])
            end
        case -2
            error('Master problem was found to be infeasible')
        case -3
            error('Master problem was found to be unbounded')
        case 6
            disp('Inaccurate solution found.  Do not update zlower')
        otherwise
            if exist('cvx_status','var')
                error(['cvx_statis is ' cvx_status ', error code is ' num2str(exitflag)])
            else
                error(['Unknown error code' num2str(exitflag)])
            end
    end
    
%     Update the upper bound on z and get the next cuts
    if opt_obj([x0;lambda0;mu0;theta0],c,N,Nbar) < zupper && ...
            newSolution
        zupper = opt_obj([x0;lambda0;mu0;theta0],c,N,Nbar);
        disp(['New best solution, updating zupper = ' num2str(zupper)])
    end
    pWorst = lambda0*numscen./(mu0-indivScens');
    disp(['Total probability = ' num2str(sum(pWorst))])
    if notPrimalSolnFound
        notPrimalSolnFound = zupper - zlower >= tolerance*min(abs(zupper),abs(zlower));
        disp(['notPrimalSolnFound = ' num2str(notPrimalSolnFound)])
    else
        disp('Primal tolerances already reached')
    end
    
    disp(' ')
end
pmle = numscen./N;

% Return costs to original scaling
c = c/scale;
indivScens = indivScens/scale;
lambda0 = lambda0/scale;
mu0 = mu0/scale;
zlower = zlower/scale;
zupper = zupper/scale;

disp(['Time elapsed = ' num2str(toc)])
read_results(x0,c,periods1)
disp(['lambda = ' num2str(lambda0) ', mu = ' num2str(mu0)])
disp(['First-stage cost = ' num2str(c*x0)])
disp(['Scenario costs = ' num2str(indivScens')])
disp(['Worst-case probabilities = ' num2str(pWorst)])
disp(['Total Probability = ' num2str(sum(pWorst))])
disp(['Gamma prime = ' num2str(gammaprime)])
disp(['Worst-case relative likelihood = ' ...
    num2str( exp(sum(numscen.*(log(pWorst)-log(pmle)))) )])
disp(['Worst-case corrected relative likelihood = ' ...
    num2str( exp(sum(numscen.*(log(pWorst./sum(pWorst))-log(pmle)))) )])
disp(['zlower = ' num2str(zlower) ', zupper = ' num2str(zupper)])
disp(['Relative error = ' num2str((zupper - zlower)/min(abs(zupper),abs(zlower)))])
disp(['Tolerance = ' num2str(tolerance)])
if zlower > zupper
    error('zlower > zupper')
end

% ------------------------------------------------------------------------
% ---------------- Accessory Functions -----------------------------------
% ------------------------------------------------------------------------

function [objective] = linear_obj(c,Nbar)
objective = [c Nbar 1 1];

% Objective function for the optimization problem
function [obj deriv] = opt_obj(x,c,N,Nbar)
obj = linear_obj(c,Nbar)*x;
deriv = linear_obj(c,Nbar);

% Get the objective value
function obj = get_first_stage_obj(x,lambda,mu,c,N,Nbar)
obj = linear_obj(c,Nbar)*[x;lambda;mu;0];
% obj = zeros(size(x));
% for ii=1:length(x)
%     obj(ii) = opt_obj([x(ii);lambda;mu;0],c,N,Nbar);
% end

% Outputs an array of h_i'(x), slope and intercept for every scenario
function [objs slope intercept] = solve_scens(x, numscen)
objs = zeros(length(numscen),1);
slope = zeros(length(numscen),length(x));
intercept = zeros(length(numscen),1);
for ii = 1:length(numscen)
    [objs(ii) slope(ii,:) intercept(ii)] = h(x,ii);
end

% Outputs a plot of the objective function in terms of x
function y = get_obj(x,lambda,mu,c,numscen,N,Nbar)
y = zeros(length(numscen),1);
for ii = 1:length(numscen)
    y(ii) = h(x,ii);
end
y = get_first_stage_obj(x,lambda,mu,c,N,Nbar) + get_exp_h(y,lambda,mu,numscen,N);
% y = get_exp_h(y,lambda,mu,numscen,N);
% y = c*x + Nbar*lambda + N*lambda*log(lambda) + mu + (numscen/N)*(-N*lambda*log(mu -y)) ;

% Find the next cut
function [expy slope intercept y] = get_cut(x,lambda,mu,numscen,N,Nbar,y,slope,intercept)
% [y slope intercept] = solve_scens(x,numscen);
intermediateSlope = zeros(size(slope) + [0,2]);
intermediateIntercept = zeros(size(intercept));
for ii=1:length(numscen)
    % October 10, 2012:
    % Note: the code below differs from what I wrote up in my preparation
    % of LRSLP-2 because I have folded the N*lambda*log(lambda) term into
    % the second stage.  Thus the slope of lambda and the intercept need to
    % be corrected.  Also, the term (y(ii) - slope(ii,:)*x) in the original
    % intermediateIntercept is wrong.  I have corrected it to be 
    % (mu - slope(ii,:)*x).  The original versions are commented out.
%     intermediateSlope(ii,:) = [(N*lambda/(mu - y(ii)))*slope(ii,:), ...
%         N*log(lambda) + N - N*log(mu - y(ii)), ...
%         -N*lambda/(mu-y(ii))];
%     intermediateIntercept(ii) = N*lambda/(mu-y(ii))*(y(ii)-slope(ii,:)*x);
    % intermediateSlope(ii,:) = [(N*lambda/(mu - y(ii)))*slope(ii,:), ...
    intermediateSlope(ii,:) = [N*lambda*(slope(ii,:)./(mu - y(ii))), ...
        N*log(lambda) + N - N*log(mu - y(ii)), ...
        -N*lambda/(mu-y(ii))];
%     intermediateIntercept(ii) = N*lambda/(mu-y(ii))*(y(ii)-slope(ii,:)*x);
end
expy = get_exp_h(y,lambda,mu,numscen,N);
% slope = [numscen/N*slope, 0, 0];
% intercept = numscen/N*intercept;
slope = numscen/N*intermediateSlope;
% intercept = numscen/N*intermediateIntercept;
intercept = expy - slope*[x;lambda;mu];

% Returns the true expected value of h for list of second stage costs y
function eh = get_exp_h(y,lambda,mu,numscen,N)
eh = numscen/N*(N*lambda*log(lambda) - N*lambda*log(mu-y));

function [feasSlope feasInt] = get_feas_cut(slopeIn,interceptIn,hIndex)
feasSlope = [slopeIn(hIndex,:), 0, -1, 0];
feasInt = interceptIn(hIndex,:);
