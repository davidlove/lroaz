function TestTwoStageGraphical

simpleLP = InitializeSimpleTwoStageLP();

gp = 0.5;
obs = ceil(50*rand(1,4));

lrlp = LRLP( simpleLP, gp, obs );

while lrlp.currentObjectiveTolerance > lrlp.objectiveTolerance
    lrlp.SolveMasterProblem;
    lrlp.SolveSubProblems;
    lrlp.GenerateCuts;
    lrlp.Plot;
%     lrlp.UpdateTrustRegionSize;
    lrlp.UpdateSolutions;
    lrlp.UpdateTolerances;
    lrlp.WriteProgress;
    pause(0.5)
end