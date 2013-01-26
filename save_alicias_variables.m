function save_alicias_variables(folder)

if nargin < 1
    folder = {'all_scenarios_new/9'};
end
startDir = pwd;
for ii=1:length(folder)
    cd(folder{ii});
    
    RunAndSave();
    
    cd(startDir)
end

function RunAndSave()

LP_Optimization;
save('orig_variables.mat','A','A_st','A_lag', ...
    'A_full','Final', ...
    'b_vec','LB','UB','Cost','Q','fval', ...
    'Var_name','Nc','Nr');