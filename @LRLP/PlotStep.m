%     Plot, to see the algorithm working
function PlotStep( obj, inVariableNumber )
x0 = obj.GetDecisions( obj.candidateSolution );
x0 = x0( inVariableNumber );
lambda0 = obj.candidateSolution.Lambda;
mu0 = obj.candidateSolution.Mu;

cMaster = obj.GetMasterc();
lMaster = obj.GetMasterl();
uMaster = obj.GetMasteru();


numPlottedPoints = 41;
switch inVariableNumber
    case obj.LAMBDA
        varPlot = linspace(0,max(0.2,2*lambda0),numPlottedPoints);
        var = '\lambda';
    case obj.MU
%         varPlot = linspace(0,max(10,2*mu0),numPlottedPoints);
        varPlot = linspace( lMaster(inVariableNumber)-1, ...
            uMaster(inVariableNumber)+1, ...
            numPlottedPoints );
        var = '\mu';
    case obj.THETA
        error('Do not plot with theta')
    otherwise
        if lMaster(inVariableNumber) == -Inf || uMaster(inVariableNumber) == Inf
            error('Plot a variable with finite bounds')
        end
        lb = obj.lpModel.l(inVariableNumber);
        ub = obj.lpModel.u(inVariableNumber);
        varPlot = linspace( max( lMaster(inVariableNumber)-1, lb ), ...
            min( uMaster(inVariableNumber)+1, ub ), ...
            numPlottedPoints );
        var = ['x[' num2str(inVariableNumber) ']'];
end

yPlot = zeros(size(varPlot));
% origCMIF = obj.candidateMuIsFeasible;

plotSolution = obj.candidateSolution.copy;
solnVector = obj.GetDecisions(plotSolution);
origSecondStageValues = obj.candidateSolution.SecondStageValues;

% Set to feasible for plotting lambda and mu
% obj.candidateMuIsFeasible = true;
for ii=1:length(varPlot)
    solnVector(inVariableNumber) = varPlot(ii);
    plotSolution.Reset;
    plotSolution.SetX( solnVector( 1:end-2-length(obj.THETA) ) );
    plotSolution.SetLambda( solnVector( obj.LAMBDA ) );
    plotSolution.SetMu( solnVector( obj.MU ) );
    plotSolution.SetTheta( solnVector( obj.THETA ), 'master' );
    
    if inVariableNumber ~= obj.LAMBDA && inVariableNumber ~= obj.MU
        obj.SolveSubProblems( plotSolution );
    else
        % Second stage values don't change, can save on recomputing
        for ss=1:obj.lpModel.numScenarios
            plotSolution.SetSecondStageValue( ss, origSecondStageValues(ss) );
        end
    end
    
    if isempty( plotSolution.MuFeasible )
        error('Must first determine feasibility of mu')
    end
    
    if plotSolution.MuFeasible
        obj.FindExpectedSecondStage( plotSolution );
%         plotSolution(obj.THETA) = obj.thetaTrue;
        yPlot(ii) = cMaster * obj.GetDecisions( plotSolution, 'true' );
    else
        yPlot(ii) = 0+1i;
    end
    
%     yplot(ii) = get_obj(xplot(ii),lambda0,mu0,cMaster,numscen,N,Nbar);
end
% Return LRLP to original status before the plot step
% obj.SolveSubProblems();
% obj.candidateMuIsFeasible = origCMIF;
% obj.FindExpectedSecondStage();

figure( inVariableNumber )
plot(varPlot(imag(yPlot)==0),yPlot(imag(yPlot)==0), ...
    varPlot(imag(yPlot)~=0),abs(imag(yPlot(imag(yPlot)~=0))),'b.', 'LineWidth',2)

yl = ylim;
% if obj.zLower > -Inf
%     yl(1) = (obj.zLower-0.1*yl(2))/(1-0.1);
    yl(1) = (cMaster*obj.GetDecisions(obj.candidateSolution)-0.1*yl(2))/(1-0.1);
% end
% if zupper < Inf && zlower > -Inf
%     a = 0.85;
%     b = 0.075;
%     zhigh = max(zupper,zlower);
%     zlow = min(zupper,zlower);
%     yl = [zhigh + (a-1)*(zhigh-zlow)/b, zhigh + a*(zhigh-zlow)/b];
% end
hold on;

if length(obj.THETA) == 1
    for cc = 1:obj.NumObjectiveCuts()
        obj.PlotObjectiveCut( inVariableNumber, cc, varPlot );
    end
end
obj.PlotBoundFunction( inVariableNumber, varPlot )

for cc=1:obj.NumFeasibilityCuts()
    obj.PlotFeasibilityCut( inVariableNumber, cc, yl );
end

% If lower bound hasn't been updated yet
lowerBound = cMaster*obj.GetDecisions(obj.candidateSolution,'master');
plot(x0, lowerBound,'k*', 'MarkerSize',8)
plot(varPlot([1,end]),[obj.zUpper obj.zUpper],'k', ...
    [x0 x0],[lowerBound obj.zUpper],'r', 'LineWidth',2)

% Plot trust region
plot(obj.trustRegionLower(inVariableNumber)*[1,1], yl, 'k--', ...
    obj.trustRegionUpper(inVariableNumber)*[1,1], yl, 'k--', ...
    'LineWidth',1)

xlim(varPlot([1,end]))
ylim(yl)
xlabel(['Decision Variable ' var], 'FontSize',14)
ylabel('Objective Function & Cuts', 'FontSize',14)
title(['\lambda = ' num2str(lambda0) ', \mu = ' num2str(mu0)], 'FontSize',14)

hold off
end
