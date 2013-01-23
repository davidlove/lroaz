if ~exist('check_folder','var')
    error('Variable ''check_folder'' must be defined');
end

if check_folder(end) ~= '/'
    check_folder = [check_folder '/'];
end
Originals = load([check_folder 'orig_variables.mat']);

disp(['Difference in A = '     num2str( max(max(abs(A-Originals.A         ))) )])
disp(['Difference in A_st = '  num2str( max(max(abs(A_st-Originals.A_st   ))) )])
disp(['Difference in A_lag = ' num2str( max(max(abs(A_lag-Originals.A_lag ))) )])
disp(['Difference in b_vec = ' num2str( max(max(abs(b_vec-Originals.b_vec ))) )])
disp(['Difference in LB = '    num2str( max(max(abs(LB-Originals.LB       ))) )])
disp(['Difference in UB = '    num2str( max(max(abs(UB-Originals.UB       ))) )])
disp(['Difference in Cost = '  num2str( max(max(abs(Cost-Originals.Cost   ))) )])
disp(['Difference in Q = '     num2str( max(max(abs(Q-Originals.Q         ))) )])
disp(['Difference in fval = '  num2str( max(max(abs(fval-Originals.fval   ))) )])

disp(['Difference in Nc = '       num2str( sum(1-strcmp(Nc,Originals.Nc       )) )])
disp(['Difference in Nr = '       num2str( sum(1-strcmp(Nr,Originals.Nr       )) )])
disp(['Difference in Var_name = ' num2str( sum(1-strcmp(Var_name,Originals.Var_name       )) )])