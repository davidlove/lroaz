function TestTwoStageGraphical

simpleLP = InitializeSimpleTwoStageLP();

gp = 0.5;
obs = ceil(50*rand(1,4));

lrlp = LRLP( simpleLP, gp, obs );

while lrlp.currentObjectiveTolerance > lrlp.objectiveTolerance
    exitFlag = lrlp.SolveMasterProblem;
    if exitFlag ~= 1
        disp(['exitFlag = ' num2str(exitFlag)])
    end
    lrlp.SolveSubProblems;
    lrlp.GenerateCuts;
%     lrlp.UpdateTrustRegionSize;
    lrlp.UpdateSolutions;
    lrlp.UpdateTolerances;
    lrlp.WriteProgress;
%     lrlp.Plot;
    pause(0.5)
end