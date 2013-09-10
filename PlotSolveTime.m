function data = PlotSolveTime( phiType, lpModel, scens )
%PlotSolveTime Plot time required to solve a phi-divergence problem vs
%number of scenarios

phi = PhiDivergence( phiType );

data = struct;

for ii = 1:length(scens)
    lp = PruneScenarios( lpModel, 1:scens(ii) );
    obs = ones(1,scens(ii));
    
    timeStart = tic;
    [solvedLRLP,c1,n1] = SolveLRLP( lp, phi, obs, -1 );
    timeIndiv = toc(timeStart);
    
    x1 = solvedLRLP.bestSolution.X;
    m1 = solvedLRLP.bestSolution.Mu;
    l1 = solvedLRLP.bestSolution.Lambda;
    r1 = solvedLRLP.calculatedDivergence;
    z1 = solvedLRLP.ObjectiveValue;
    
    if ii==1
        data.scens = scens;
        data.x = zeros(length(x1),length(scens));
        data.lambda = zeros(length(l1),length(scens));
        data.mu = zeros(length(m1),length(scens));
        data.objVals = zeros(length(z1),length(scens));
        data.calcRho = zeros(length(r1),length(scens));
        data.numProbs = zeros(1,length(scens));
        data.numCuts = zeros(1,length(scens));
        data.timeRuns = zeros(1,length(scens));
    end
    
    data.x(:,ii) = x1;
    data.lambda(:,ii) = l1;
    data.mu(:,ii) = m1;
    data.objVals(:,ii) = z1;
    data.calcRho(:,ii) = r1;
    data.numProbs(ii) = n1;
    data.numCuts(ii) = c1;
    data.timeRuns(ii) = timeIndiv;
end

plot(data.scens, data.timeRuns, 'o', 'MarkerSize',10)
xlabel( 'N', 'FontSize',16 )
ylabel( 'Solution Time', 'FontSize',16)

end

