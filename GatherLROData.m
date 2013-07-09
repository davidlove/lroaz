function GatherLROData( phiType, saveFileName )

if nargin < 2
    saveFileName = 'saved_variables.mat';
    if nargin < 1
        phiType = 'kl';
    end
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
phi = PhiDivergence( phiType );

obs = 1*[1 1 1 1];

if isinf(phi.limit) && any(obs == 0)
    error(['Phi divergence ' phi.divergence ' does not allow poping up probabilities'])
end

alpha = 10.^linspace(-2.93,0,50);
alpha(end) = mean(alpha(end-1:end));
rho = phi.SecondDerivativeAt1() / (2*sum(obs)) * ...
    chi2inv(1-alpha,lpModel.numScenarios - 1);

iiSet = 1:length(rho);

for ii = iiSet
    timeStart = tic;
    [solvedLRLP,c1,n1] = SolveLRLP( lpModel, phi, obs, rho(ii), 'multi' );
    timeIndiv = toc(timeStart);
    
    p1 = solvedLRLP.pWorst;
    x1 = solvedLRLP.bestSolution.X;
    m1 = solvedLRLP.bestSolution.Mu;
    l1 = solvedLRLP.bestSolution.Lambda;
    s1 = solvedLRLP.bestSolution.SecondStageValues;
    r1 = solvedLRLP.calculatedDivergence;
    z1 = solvedLRLP.ObjectiveValue;
    
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
    
    save(saveFileName,'rho','pWorst','x','lambda','mu', ...
        'scenCosts','objVals', ...
        'calcRho', 'alpha', ...
        'numProbs','numCuts','timeRuns');
%         'exitFlags')

    clear solvedLRLP
end

