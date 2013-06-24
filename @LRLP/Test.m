function Test( obj )
%TEST Run tests on LRLP class

% Ensure that SolveMasterProblem runs -- at one point, phi divergences with
% infinite limit would not initialize Mu correctly and not deal with the
% feasibility matrix correctly.
% NOTE: SolveMasterProblem must be called BEFORE I start changing values of
% x, lambda, and mu.
obj.SolveMasterProblem

% lp parameters
n = obj.lpModel.numScenarios;

% phi parameters
divergence = obj.phi.divergence;
phi2deriv = obj.phi.SecondDerivativeAt1();
limit = obj.phi.limit;

% LRLP Parameters
N = obj.numObsTotal;
q = obj.numObsPerScen/N;

% IMPORTANT: These tests only work for default value of rho, with multicut!

[c] = obj.GetMasterc();
assertElementsAlmostEqual( c(1:length(obj.lpModel.c)), obj.lpModel.c )
assertElementsAlmostEqual( c(obj.LAMBDA), phi2deriv/(2*N)*chi2inv(0.95,n-1) )
assertElementsAlmostEqual( c(obj.MU), 1 )
assertElementsAlmostEqual( c(obj.THETA), q )

% Reset parameters set by InitializeBenders
obj.candidateSolution.Reset()
obj.objectiveCutsMatrix = [];
obj.objectiveCutsRHS = [];
obj.feasibilityCutsMatrix = [];
obj.feasibilityCutsRHS = [];

% Select values for phi*, choose second stage values to satisfy these s
% values
x = zeros(1,length(obj.lpModel.c));
lambda = 3;
mu = 1/2;
switch divergence
    case 'burg'
        conjvals = linspace(-3,0.5,n)';
        s = 1 - exp(-conjvals);
    case 'kl'
        conjvals = linspace(-0.5,5,n)';
        s = log(conjvals+1);
    case 'chi2'
        conjvals = linspace(-3,1.5,n)';
        s = 1 - (1 - conjvals/2).^2;
    case 'hellinger'
        conjvals = linspace(-0.5,5,n)';
        s = conjvals./(1+conjvals);
    otherwise
        error(['Test not set up for divergence: ' divergence])
end
h = mu + s*lambda;

% Set up to test finding expected second stage
obj.candidateSolution.SetX(x)
obj.candidateSolution.SetLambda(lambda)
obj.candidateSolution.SetMu(mu)
for ii = 1:n
    obj.candidateSolution.SetSecondStageValue(ii,h(ii))
end
predictedThetaTrue = lambda*conjvals;
obj.FindExpectedSecondStage( obj.candidateSolution )
assertElementsAlmostEqual( obj.candidateSolution.ThetaTrue(), predictedThetaTrue )

% Select values for phi*', choose second stage values to satisfy these s
% values
switch divergence
    case 'burg'
        conjderivs = linspace(0.1,5,n)';
        s = 1 - 1./conjderivs;
    case 'kl'
        conjderivs = linspace(0.1,5,n)';
        s = log(conjderivs);
    case 'chi2'
        conjderivs = linspace(0.1,5,n)';
        s = 1 - 1./(conjderivs.^2);
    case 'hellinger'
        conjderivs = linspace(0.1,5,n)';
        s = 1 - 1./sqrt(conjderivs);
    otherwise
        error(['Test not set up for divergence: ' divergence])
end
h = mu + s*lambda;
for ii = 1:n
    obj.candidateSolution.SetSecondStageValue(ii,h(ii))
end
conjvals = obj.phi.Conjugate(s);

% Define second stage dual values
secondStageSlope = repmat( (1:n)', 1, length(x) );
secondStageInt = -(1:n)';
for ii=1:n
    obj.candidateSolution.SetSecondStageDual(ii, secondStageSlope(ii), 'slope')
    obj.candidateSolution.SetSecondStageDual(ii, secondStageInt(ii), 'int')
end
ident = eye(n);
predictedCutsMatrix = zeros(n,length(c));
for ii=1:n
    predictedCutsMatrix(ii,:) = [conjderivs(ii)* secondStageSlope(ii), ...
        conjvals(ii) - conjderivs(ii)*s(ii), ...
        -conjderivs(ii), ...
        -ident(ii,:)];
end
predictedCutsMatrix = sparse(predictedCutsMatrix);
predictedCutsRHS = predictedCutsMatrix*[x;lambda;mu;zeros(n,1)] - predictedThetaTrue;

obj.GenerateObjectiveCut()

assertElementsAlmostEqual( obj.objectiveCutsMatrix, predictedCutsMatrix )
assertElementsAlmostEqual( obj.objectiveCutsRHS, predictedCutsRHS )

% Make variables infeasible, test feasibility cut, find feasible mu
if limit < Inf
    f = 1;
    s(f) = limit + 1;
    h = mu + s*lambda;
    for ii = 1:n
        obj.candidateSolution.SetSecondStageValue(ii,h(ii))
    end
    predictedFeasMatrix = sparse([secondStageSlope(f), -limit, -1, zeros(1,n)]);
    predictedFeasRHS = -secondStageInt(f);
    
    obj.GenerateFeasibilityCut()
    assertElementsAlmostEqual( obj.feasibilityCutsMatrix, predictedFeasMatrix )
    assertElementsAlmostEqual( obj.feasibilityCutsRHS, predictedFeasRHS )
    
    obj.FindFeasibleMu()
    assertTrue( all( obj.candidateSolution.S() < limit ) )
end

obj.CalculateProbability

end