function GatherLROData(saveFileName )

if nargin < 1
    saveFileName = 'saved_variables.mat';
end

% lpModel = InitializeSimpleTwoStageLP();

waterFolders = { ...
                 'SP C/5/', ...
                 'SP C/6/', ...
                 'SP C/7/', ...
                 'SP C/8/', ...
               };
years = 41;
timeLag = 1;
firstStageYears = 5;

lpModel = LPModel( waterFolders, years, timeLag, firstStageYears );

obs = 1*[1 1 1 1];

dgp = 0.01;
gp = dgp:dgp:1-dgp;

iiSet = 1:length(gp);

for ii = iiSet
    timeStart = tic;
    [solvedLRLP,c1,n1] = SolveLRLP( lpModel, gp(ii), obs );
    timeIndiv = toc(timeStart);
    
    p1 = solvedLRLP.pWorst;
    x1 = solvedLRLP.bestSolution.X;
    m1 = solvedLRLP.bestSolution.Mu;
    l1 = solvedLRLP.bestSolution.Lambda;
    s1 = solvedLRLP.bestSolution.SecondStageValues;
    r1 = solvedLRLP.relativeLikelihood;
    
    if ii==1
        pWorst = zeros(length(p1),length(gp));
        x = zeros(length(x1),length(gp));
        lambda = zeros(length(l1),length(gp));
        mu = zeros(length(m1),length(gp));
        scenCosts = zeros(length(s1),length(gp));
%         zCost = zeros(length(z1),length(gp));
        likelihood = zeros(length(r1),length(gp));
%         exitFlags = zeros(1,length(gp));
        numProbs = zeros(1,length(gp));
        numCuts = zeros(1,length(gp));
        timeRuns = zeros(1,length(gp));
    end
    
    pWorst(:,ii) = p1';
    x(:,ii) = x1;
    lambda(:,ii) = l1;
    mu(:,ii) = m1;
    scenCosts(:,ii) =s1;
%     zCost(:,ii) = z1;
    likelihood(:,ii) = r1;
%     exitFlags(ii) = eF;
    numProbs(ii) = n1;
    numCuts(ii) = c1;
    timeRuns(ii) = timeIndiv;
    
    save(saveFileName,'gp','pWorst','x','lambda','mu', ...
        'scenCosts',...%'zCost',
        'likelihood', ...
        'numProbs','numCuts','timeRuns');
%         'exitFlags')

    clear solvedLRLP
end

