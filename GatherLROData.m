function GatherLROData( varargin )
% GatherLROData solves the phi-divergence model for multiple values of rho
% and reports the results
%
% Keyword arguments (*required arguments):
%    phi*: PhiDivergence object or string defining the phi-divergence
%    alpha: Asymptotical confidence levels.  Conflicts with rho
%    lp: LPModel object defining the SLP-2
%    obs: vector of observations numbers of the scenarios
%    restart: true if restarting from previously saved results
%    rho: Confidence region size rho.  Conflicts with alpha
%    savefile: String or file object denoting the file to save data to

% -------------------------------------------------------------------
% Input Parsing
% -------------------------------------------------------------------

requiredArgs = {'phi'};

artag = false;

if mod(length(varargin),2) == 1
    error('Arguments must be key, value pairs')
end

for vv = 1:2:length(varargin)
    key = varargin{vv};
    value = varargin{vv+1};
    switch key
        case 'alpha'
            if artag
                error('Only accepts one of alpha, rho')
            elseif exist('restart', 'var') && restart
                warning('Restarting ignores alpha and rho')
            end
            artag = true;
            alpha = value;
        case 'lp'
            lp = value;
        case 'obs'
            obs = value;
        case 'phi'
            if isa(value, 'PhiDivergence')
                phi = value;
            elseif ischar(value)
                phi = PhiDivergence( value );
            else
                error('Phi must be PhiDivergence object or string')
            end
        case 'restart'
            restart = value;
        case 'rho'
            if artag
                error('Only accepts one of alpha, rho')
            elseif exist('restart', 'var') && restart
                warning('Restarting ignores alpha and rho')
            end
            artag = true;
            rho = value;
        case 'savefile'
            saveFileName = value;
        otherwise
            error(['Unknown variable ', key])
    end
end

% -------------------------------------------------------------------
% Default variable assignments
% -------------------------------------------------------------------

if ~exist('lp', 'var')
        waterFolders = { ...
            'SP C/5/', ...
            'SP C/6/', ...
            'SP C/7/', ...
            'SP C/8/', ...
            };
        years = 41;
        timeLag = 1;
        firstStageYears = 5;
        
        lp = LPModel( waterFolders, years, timeLag, firstStageYears );
end
if ~exist('obs', 'var')
    obs = ones(1,lp.numScenarios);
end
if ~exist('restart', 'var')
    restart = false;
end
if ~exist('saveFileName', 'var')
    saveFileName = 'saved_variables.mat';
end

if ~restart
    if ~exist('rho', 'var')
        if ~exist('alpha', 'var')
            alpha = 10.^linspace(-2.93,0,50);
            alpha(end) = mean(alpha(end-1:end));
        end
        rho = phi.Rho(alpha, obs);
    elseif ~exist('alpha', 'var')
        alpha = phi.Alpha(rho, obs);
    end
    
    iiSet = 1:length(rho);
else
    load( saveFileName );
    
    iiSet = find( max(pWorst) == 0 );
end

% -------------------------------------------------------------------
% Required Variables
% -------------------------------------------------------------------

for aa = requiredArgs
    if ~exist(aa{1}, 'var')
        error([aa{1}, ' is required but not defined'])
    end
end

% -------------------------------------------------------------------
% Value Checking
% -------------------------------------------------------------------

if ~isa(lp, 'LPModel')
    error('lp must be an LPModel object')
end

if isinf(phi.limit) && any(obs == 0)
    error(['Phi divergence ' phi.divergence ' does not allow poping up probabilities'])
end
if restart ~= true && restart ~= false
    error('restart must be a boolean value')
end

% -------------------------------------------------------------------
% -------------------------------------------------------------------



for ii = iiSet
    timeStart = tic;
    
    try
        [solvedLRLP,c1,n1] = SolveLRLP( 'lp',lp, 'phi',phi, 'obs',obs, ...
            'rho',rho(ii), 'cuttype','multi' );
    catch err
        warning('Error detected in SolveLRLP')
        getReport(err)
        pause(30)
        continue
    end
    
    timeIndiv = toc(timeStart);
    
    p1 = solvedLRLP.pWorst;
    x1 = solvedLRLP.bestSolution.X;
    m1 = solvedLRLP.bestSolution.Mu;
    l1 = solvedLRLP.bestSolution.Lambda;
    s1 = solvedLRLP.bestSolution.SecondStageValues / solvedLRLP.objectiveScale;
    r1 = solvedLRLP.calculatedDivergence;
    z1 = solvedLRLP.ObjectiveValue / solvedLRLP.objectiveScale;
    y1 = cell(lp.numScenarios, 1);
    for ss = 1:length(y1)
        y1{ss} = solvedLRLP.bestSolution.SecondStageSolution(ss);
    end
    
    if ii==1
        pWorst = zeros(length(p1),length(rho));
        x = zeros(length(x1),length(rho));
        lambda = zeros(length(l1),length(rho));
        mu = zeros(length(m1),length(rho));
        scenCosts = zeros(length(s1),length(rho));
        objVals = zeros(length(z1),length(rho));
        calcRho = zeros(length(r1),length(rho));
%         exitFlags = zeros(1,length(gp));
        numProbs = zeros(1,length(rho));
        numCuts = zeros(1,length(rho));
        timeRuns = zeros(1,length(rho));
        scenSolns = cell(length(y1),length(rho));
        profileNames = lp.profileNames;
    end
    
    pWorst(:,ii) = p1';
    x(:,ii) = x1;
    lambda(:,ii) = l1;
    mu(:,ii) = m1;
    scenCosts(:,ii) =s1;
    objVals(:,ii) = z1;
    calcRho(:,ii) = r1;
%     exitFlags(ii) = eF;
    numProbs(ii) = n1;
    numCuts(ii) = c1;
    timeRuns(ii) = timeIndiv;
    scenSolns(:,ii) = y1;
    
    save(saveFileName,'rho','pWorst','x','lambda','mu', ...
        'scenCosts','scenSolns','objVals', ...
        'calcRho', 'alpha', ...
        'numProbs','numCuts','timeRuns', ...
        'profileNames');
%         'exitFlags')

    clear solvedLRLP
end

