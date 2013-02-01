function [lrlp, outTotalCuts, outTotalProbs] = ...
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
            case {-2,-3,-4,-5}
                % Do nothing extra
            otherwise
                error( ['Unknown error code: ' num2str(exitFlag)] )
        end
        lrlp.DeleteOldestCut;
        if lrlp.NumFeasibilityCuts > 1
            lrlp.DeleteOldestFeasibilityCut;
        end
        disp([num2str(lrlp.NumObjectiveCuts) ' Objective Cuts Remaining, ' ...
            num2str(lrlp.NumFeasibilityCuts) ' Feasibility Cuts Remaining.'])
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

outTotalCuts = totalCutsMade;
outTotalProbs = totalProblemsSolved;

end