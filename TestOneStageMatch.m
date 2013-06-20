function TestOneStageMatch

dirs = 1:16;
Period = 41;
timeLag = 1;

for ii=1:length(dirs)
    inputLocation = strcat('SP C/',num2str(dirs(ii)),'/');
    clear lp;
    lp = LPModel(inputLocation,Period,timeLag);
    [Q,fval] = linprog( lp.c, [], [], lp.A, lp.b, lp.l, lp.u );
%     Q = roundn(Q',-2);
    Q = lp.ReadResults(Q,inputLocation);
    
    disp(' ')
    disp(['Checking folder: ' inputLocation])
    disp(' ')
    alicia = load([inputLocation,'orig_variables.mat']);
    % assertEqual does not work when comparing sparse (lp.A) with dense
    % (alicia.A_full) matrices.
    assert(isequal(lp.A,alicia.A_full))
    assertEqual(lp.b,alicia.b_vec)
    assertEqual(lp.c,alicia.Cost')
    assertEqual(lp.l,alicia.LB)
    assertEqual(lp.u,alicia.UB)
    assertEqual(Q,alicia.Q)
    assertEqual(fval,alicia.fval)
    assertEqual(lp.Final,alicia.Final)
    
    disp(['Folder ' inputLocation ' passed tests'])
    disp(' ')
end

