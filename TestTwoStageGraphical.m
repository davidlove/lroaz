function TestTwoStageGraphical

simpleLP = InitializeSimpleTwoStageLP();
phi = PhiDivergence('lro');
origLP = simpleLP;
assert( isequal( simpleLP, origLP ) )

problemCases = [];

numKnownProblems = size(problemCases,1);
numNewProblems = 10;
numProblems = numKnownProblems + numNewProblems;

cuts = zeros(numProblems,1);
probs = cuts;
tols = cuts;
ptols = cuts;

for ii = 1:numProblems
    if ii <= numKnownProblems
        rho = problemCases(ii,1);
        obs = problemCases(ii,2:end);
    else
        obs = ceil(50*rand(1,4));
        rho = -1;
    end
    
    [solvedLRLP,cuts(ii),probs(ii)] = SolveLRLP( simpleLP, phi, obs, rho, 'multi', false );
    
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
disp('   Maximums')
disp([max(tols),max(ptols)])
end