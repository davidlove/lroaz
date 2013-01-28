function TestOneStageTimeLag

Period = 2;
timeLag = 12;
inputLocation = 'all_scenarios_hwee/';
clear lp;
lp = LPModel(inputLocation,Period,timeLag);
[Q fval] = linprog( lp.c, [], [], lp.A, lp.b, lp.l, lp.u );
%     Q = roundn(Q',-2);
Q = lp.ReadResults(Q,inputLocation);

hwee = load([inputLocation,'orig_variables.mat']);
assertEqual(lp.b,hwee.b_vec)
assertEqual(lp.c,hwee.Cost')
assertEqual(lp.l,hwee.LB)
assertEqual(lp.u,hwee.UB)
% assert(nnz(lp.A ~= hwee.A_full) == 0)
% assertEqual(Q,hwee.Q)
% assertEqual(fval,hwee.fval)
% assertEqual(lp.Final,hwee.Final)

