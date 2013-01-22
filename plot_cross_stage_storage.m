function plot_cross_stage_storage ...
    ( solutionsFile,connectionsFile,cellInputFile )

% plot_cross_stage_storage plots how the amount of water stored at the end
% of the first stage varies with gamma'

clear get_stage_vectors

if nargin < 3
    cellInputFile = { ...
        'all_scenarios/8/Inputs.xlsx', ...
        'all_scenarios/5/Inputs.xlsx', ...
        'all_scenarios/6/Inputs.xlsx', ...
        'all_scenarios/7/Inputs.xlsx', ...
        };
    if nargin < 2
        connectionsFile = 'all_scenarios/8/Connections.xlsx';
        if nargin < 1
            solutionsFile = 'saved_variables_01_99.mat';
        end
    end
end

solns = load(solutionsFile);
[~,~,~,~,~,B] = get_stage_vectors(2,1, ...
    connectionsFile,cellInputFile,41,10);
[~,Nc,Nr] = get_stage_vectors(-1);

good = find(solns.exitFlags==1);
gp = solns.gp;
x = solns.x;

[r,c] = find(B);

[~,cIndex] = unique(c,'first');
cUnique=c(sort(cIndex));
carryX = x(cUnique,good);
names = Nc(mod(cUnique,length(Nc)));
plot(gp(good),carryX,'.')