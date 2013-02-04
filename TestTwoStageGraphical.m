function TestTwoStageGraphical

simpleLP = InitializeSimpleTwoStageLP();
origLP = simpleLP;
assert( isequal( simpleLP, origLP ) )

% gp = 0.5;
problemCases = [0.81, 46   7  46  32; ...
                0.89,  7  20  47  46; ...
                0.63,  2  46  41  38; ...
                0.50, 35  28  22  33];
% obs = ceil(50*rand(1,4));

numKnownProblems = size(problemCases,1);
numNewProblems = 10;
numProblems = numKnownProblems + numNewProblems;

cuts = zeros(numProblems,1);
probs = cuts;
tols = cuts;
ptols = cuts;

for ii = 1:numProblems
    if ii <= numKnownProblems
        gp = problemCases(ii,1);
        obs = problemCases(ii,2:end);
    else
        gp = roundn( rand(), -2 );
        obs = ceil(50*rand(1,4));
    end
    
    [solvedLRLP,cuts(ii),probs(ii)] = SolveLRLP( simpleLP, gp, obs, false );
    
    tols(ii) = solvedLRLP.currentObjectiveTolerance;
    ptols(ii) = solvedLRLP.currentProbabilityTolerance;
    
    assert( isequal( simpleLP, origLP ) )
    
    clear solvedLRLP    
end

disp([cuts,probs;...
      max(cuts),max(probs)])
  format short e
disp([tols,ptols;...
      max(tols),max(ptols)])
end