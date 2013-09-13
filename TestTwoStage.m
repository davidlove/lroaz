function TestTwoStage

simpleLP = InitializeSimpleTwoStageLP();
phiTypes = {'burg','kl','chi2','mchi2'};
% phi = PhiDivergence('mchi2');
origLP = simpleLP;
assert( isequal( simpleLP, origLP ) )

problemCases = {'mchi2', [-1, 27, 39, 47, 7]}; % mchi2, poorly scaled lp

numKnownProblems = size(problemCases,1);
numNewProblems = 2;
numProblems = numKnownProblems + numNewProblems*length(phiTypes);

cuts = zeros(numProblems,1);
probs = cuts;
tols = cuts;
ptols = cuts;

for ii = 1:numProblems
    if ii <= numKnownProblems
        phi = PhiDivergence(problemCases{ii,1});
        details = problemCases{ii,2};
        rho = details(ii,1);
        obs = details(ii,2:end);
    else
        type = phiTypes{ceil((ii-numKnownProblems)/numProblems*length(phiTypes))};
        phi = PhiDivergence( type );
        obs = ceil(50*rand(1,4));
        rho = -1;
    end
    
    [solvedLRLP,cuts(ii),probs(ii)] = SolveLRLP( 'lp',simpleLP, 'phi',phi, ...
        'obs',obs, 'rho',rho, 'cuttype','multi', 'isgraphical',false );
    
    tols(ii) = solvedLRLP.currentObjectiveTolerance;
    ptols(ii) = solvedLRLP.currentProbabilityTolerance;
    
    assert( isequal( simpleLP, origLP ) )
    
    clear solvedLRLP    
end

disp('   ------------------------')
disp('    Cuts  Problems')
disp([cuts,probs])
disp('    Maximums')
disp([max(cuts),max(probs)])
format short e
disp('   ------------------------')
disp('   Objective    Probability')
disp([tols,ptols])

mtol = max(tols);
mp = max(ptols);
disp('   Maximums')
disp([mtol,mp])

assert(mtol <= 1e-6)
assert(mp <= 1e-2)
end