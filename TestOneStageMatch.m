function TestOneStageMatch

dirs = 1:12;
Period = 41;
timeLag = 1;

for ii=1:length(dirs)
    inputLocation = strcat('all_scenarios/',num2str(dirs(ii)),'/');
    clear lp;
    lp = LPModel(inputLocation,Period,timeLag);
    [Q fval] = linprog( lp.c, [], [], lp.A, lp.b, lp.l, lp.u );
%     Q = roundn(Q',-2);
    Q = lp.ReadResults(Q,inputLocation);
    
    alicia = load([inputLocation,'orig_variables.mat']);
    assert(nnz(lp.A ~= alicia.A_full) == 0)
    assertEqual(lp.b,alicia.b_vec)
    assertEqual(lp.c,alicia.Cost')
    assertEqual(lp.l,alicia.LB)
    assertEqual(lp.u,alicia.UB)
    assertEqual(Q,alicia.Q)
    assertEqual(fval,alicia.fval)
    assertEqual(lp.Final,alicia.Final)
end

