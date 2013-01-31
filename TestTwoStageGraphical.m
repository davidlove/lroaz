function TestTwoStageGraphical

simpleLP = InitializeSimpleTwoStageLP();

gp = 0.5;
obs = ceil(50*rand(1,4));
obs = [35  28  22  33];

lrlp = LRLP( simpleLP, gp, obs );

while lrlp.currentObjectiveTolerance > lrlp.objectiveTolerance
    exitFlag = lrlp.SolveMasterProblem;
    if exitFlag ~= 1
        disp(['exitFlag = ' num2str(exitFlag)])
    end
    cS = lrlp.candidateSolution;
    
    lrlp.SolveSubProblems;
    assert( isequal( cS, lrlp.candidateSolution ), ...
        num2str([cS, lrlp.candidateSolution]) )
    
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
    
%     lrlp.UpdateTrustRegionSize;
%     assert( isequal( cS, lrlp.candidateSolution ), ...
%         num2str([cS, lrlp.candidateSolution]) )
    
    lrlp.UpdateSolutions;
    assert( isequal( cS, lrlp.candidateSolution ), ...
        num2str([cS, lrlp.candidateSolution]) )
    
    lrlp.UpdateTolerances;
    assert( isequal( cS, lrlp.candidateSolution ), ...
        num2str([cS, lrlp.candidateSolution]) )
    
    lrlp.WriteProgress;
    assert( isequal( cS, lrlp.candidateSolution ), ...
        num2str([cS, lrlp.candidateSolution]) )
    
%     pause(.5)
end