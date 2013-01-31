function TestTwoStageGraphical

simpleLP = InitializeSimpleTwoStageLP();
origLP = simpleLP;
assert( isequal( simpleLP, origLP ) )

% gp = 0.5;
problemCases = [0.89,  7  20  47  46; ...
                0.63,  2  46  41  38; ...
                0.50, 35  28  22  33];
% obs = ceil(50*rand(1,4));

numKnownProblems = size(problemCases,1);
numNewProblems = 20;
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
    
    [cuts(ii),probs(ii),tols(ii),ptols(ii)] ...
        = SolveLRLP( simpleLP, gp, obs );
    
    assert( isequal( simpleLP, origLP ) )
    
    
end

disp([cuts,probs;...
      max(cuts),max(probs)])
  format long e
disp([tols,ptols;...
      max(tols),max(ptols)])
end

function [cutsOut, probOut, tolOut, ptolOut] = ...
    SolveLRLP( simpleLP, gp, obs )

lrlp = LRLP( simpleLP, gp, obs );

totalProblemsSolved = 1;
totalCutsMade = 1;
while lrlp.currentObjectiveTolerance > lrlp.objectiveTolerance
    totalProblemsSolved = totalProblemsSolved + 1;
    exitFlag = lrlp.SolveMasterProblem;
    if exitFlag ~= 1
        disp(['exitFlag = ' num2str(exitFlag)])
        switch exitFlag
            case 0
                lrlp.DoubleIterations;
                lrlp.DeleteOldestCut;
            case {-2,-3,-4,-5}
                lrlp.DeleteOldestCut;
            otherwise
                error( ['Unknown error code: ' num2str(exitFlag)] )
        end
        continue
    end
    cS = lrlp.candidateSolution;
    
    lrlp.SolveSubProblems;
    assert( isequal( cS, lrlp.candidateSolution ), ...
        num2str([cS, lrlp.candidateSolution]) )
    
    totalCutsMade = totalCutsMade + 1;
    lrlp.GenerateCuts;
    assert( ~lrlp.candidateMuIsFeasible ...
        || isequal( cS, lrlp.candidateSolution ), ...
        num2str([cS, lrlp.candidateSolution]) )
    cS = lrlp.candidateSolution;
    
    % Note, putting plot after updating trust region or solution will
    % result in it plotting the wrong trust region
    %     lrlp.Plot;
    %     assert( isequal( cS, lrlp.candidateSolution ), ...
    %         num2str([cS, lrlp.candidateSolution]) )
    
    lrlp.UpdateTrustRegionSize;
    assert( isequal( cS, lrlp.candidateSolution ), ...
        num2str([cS, lrlp.candidateSolution]) )
    
    lrlp.UpdateSolutions;
    assert( isequal( cS, lrlp.candidateSolution ), ...
        num2str([cS, lrlp.candidateSolution]) )
    
    lrlp.UpdateTolerances;
    assert( isequal( cS, lrlp.candidateSolution ), ...
        num2str([cS, lrlp.candidateSolution]) )
    
    lrlp.WriteProgress;
    assert( isequal( cS, lrlp.candidateSolution ), ...
        num2str([cS, lrlp.candidateSolution]) )
    
    disp(['Total cuts made: ' num2str(totalCutsMade)])
    disp(['Total problems solved: ' num2str(totalProblemsSolved)])
    
    %     pause(.5)
end

cutsOut = totalCutsMade;
probOut = totalProblemsSolved;
tolOut = lrlp.currentObjectiveTolerance;
ptolOut = lrlp.currentProbabilityTolerance;

end