%     Plot, to see the algorithm working
function PlotStep( obj, inVariableNumber )
plotSolution = obj.candidateSolution;
x0 = plotSolution( inVariableNumber );
lambda0 = plotSolution( obj.LAMBDA );
mu0 = plotSolution( obj.MU );

cMaster = obj.GetMasterc();
lMaster = obj.GetMasterl();
uMaster = obj.GetMasteru();


numPlottedPoints = 41;
switch inVariableNumber
    case obj.LAMBDA
        varPlot = linspace(0,max(0.2,2*lambda0),numPlottedPoints);
    case obj.MU
%         varPlot = linspace(0,max(10,2*mu0),numPlottedPoints);
        varPlot = linspace( lMaster(inVariableNumber)-1, ...
            uMaster(inVariableNumber)+1, ...
            numPlottedPoints );
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
end

yPlot = zeros(size(varPlot));
for ii=1:length(varPlot)
    plotSolution(inVariableNumber) = varPlot(ii);
    
    if inVariableNumber ~= obj.LAMBDA && inVariableNumber ~= obj.MU
        obj.SolveSubProblems( plotSolution );
    end
    
    if obj.candidateMuIsFeasible
        obj.FindExpectedSecondStage( plotSolution );
        plotSolution(obj.THETA) = obj.thetaTrue;
        yPlot(ii) = cMaster * plotSolution;
    else
        yPlot(ii) = 0+1i;
    end
    
%     yplot(ii) = get_obj(xplot(ii),lambda0,mu0,cMaster,numscen,N,Nbar);
end
obj.SolveSubProblems();
obj.FindExpectedSecondStage();

figure( inVariableNumber )
plot(varPlot(imag(yPlot)==0),yPlot(imag(yPlot)==0), ...
    varPlot(imag(yPlot)~=0),abs(imag(yPlot(imag(yPlot)~=0))),'b.', 'LineWidth',2)

yl = ylim;
if obj.zLower > -Inf
    yl(1) = (obj.zLower-0.1*yl(2))/(1-0.1);
end
% if zupper < Inf && zlower > -Inf
%     a = 0.85;
%     b = 0.075;
%     zhigh = max(zupper,zlower);
%     zlow = min(zupper,zlower);
%     yl = [zhigh + (a-1)*(zhigh-zlow)/b, zhigh + a*(zhigh-zlow)/b];
% end
hold on;

for cc = 1:obj.NumObjectiveCuts()
    obj.PlotObjectiveCut( inVariableNumber, cc, varPlot );
end

for cc=1:obj.NumFeasibilityCuts()
    obj.PlotFeasibilityCut( inVariableNumber, cc, yl );
end

% plot(x0,obj.zLower,'ko', 'MarkerSize',8)
plot(x0,cMaster*obj.candidateSolution,'ko', 'MarkerSize',8)
plot(varPlot([1,end]),[obj.zUpper obj.zUpper],'k', ...
    [x0 x0],[obj.zLower obj.zUpper],'r', 'LineWidth',2)

% Plot trust region
plot(obj.trustRegionLower(inVariableNumber)*[1,1], yl, 'k:', ...
    obj.trustRegionUpper(inVariableNumber)*[1,1], yl, 'k:', ...
    'LineWidth',2)

xlim(varPlot([1,end]))
ylim(yl)
xlabel('Decision Variable x', 'FontSize',14)
ylabel('Objective Function & Cuts', 'FontSize',14)
title(['\lambda = ' num2str(lambda0) ', \mu = ' num2str(mu0)], 'FontSize',14)

hold off
end
