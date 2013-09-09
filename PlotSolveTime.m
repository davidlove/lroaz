function PlotSolveTime( phiType, lpModel, scens )
%PlotSolveTime Plot time required to solve a phi-divergence problem vs
%number of scenarios

phi = PhiDivergence( phiType );

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
        x = zeros(length(x1),length(scens));
        lambda = zeros(length(l1),length(scens));
        mu = zeros(length(m1),length(scens));
        objVals = zeros(length(z1),length(scens));
        calcRho = zeros(length(r1),length(scens));
        numProbs = zeros(1,length(scens));
        numCuts = zeros(1,length(scens));
        timeRuns = zeros(1,length(scens));
    end
    
    x(:,ii) = x1;
    lambda(:,ii) = l1;
    mu(:,ii) = m1;
    objVals(:,ii) = z1;
    calcRho(:,ii) = r1;
    numProbs(ii) = n1;
    numCuts(ii) = c1;
    timeRuns(ii) = timeIndiv;
end

plot(scens, timeRuns, 'o')

end

