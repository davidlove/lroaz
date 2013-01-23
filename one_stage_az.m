function Q = one_stage_az( inputLocation )

clear get_stage_vectors
clc

if nargin < 1
    inputLocation = 'all_scenarios_new/9/';
end

% one_stage_az uses the code from lroaz_version15 to do a one-stage problem
% that takes a given input

tic;

ConnectionsFile = [ inputLocation 'Connections.xlsx'];
cellInputFile = { [inputLocation 'Inputs.xlsx'] };
solutionFile = [inputLocation 'Solution.xlsx'];

% Problem Parameters
Period = 41;
periods1 = Period;
Time_Lag = 1;

%% Clear inputs file if connections have been altared

% Question1 = input('Has the connection matrix changed?  Yes=1  No=2 :');
% if Question1 == 1
%     % Document old variable names in excel for comparison to new
%     [~,str] = xlsread('cellInputFile{1}','b_vec','B:B');
%     str(1) = {'Old'};
%     xlswrite('cellInputFile{1}',blank,'b_vec','A');
%     xlswrite('cellInputFile{1}',str,'b_vec','A');
%     xlswrite('cellInputFile{1}',blank,'b_vec','B2');
%     
%     [~,str] = xlsread('cellInputFile{1}','Upper Bounds','B:B');
%     str(1) = {'Old'};
%     xlswrite('cellInputFile{1}',blank,'Upper Bounds','A');
%     xlswrite('cellInputFile{1}',str,'Upper Bounds','A');
%     xlswrite('cellInputFile{1}',blank,'Upper Bounds','B2');
%     
%     xlswrite('cellInputFile{1}',blank,'Lower Bounds','A');
%     xlswrite('cellInputFile{1}',str,'Lower Bounds','A');
%     xlswrite('cellInputFile{1}',blank,'Lower Bounds','B2');
%     
%     xlswrite('cellInputFile{1}',blank,'Costs','A');
%     xlswrite('cellInputFile{1}',str,'Costs','A');
%     xlswrite('cellInputFile{1}',blank,'Costs','B2');
%     
%     xlswrite('cellInputFile{1}',blank,'kWh','A');
%     xlswrite('cellInputFile{1}',str,'kWh','A');
%     xlswrite('cellInputFile{1}',blank,'kWh','B2');
%     sound(y, Fs);
% end

%% Build and solve LP

[c,A_full,b_vec,LB,UB] = get_stage_vectors(1,1, ...
    ConnectionsFile,cellInputFile,Period,periods1,Time_Lag);
[Var_name,Nc,Nr] = get_stage_vectors('names');
[A,A_st,A_lag,cost,Zones] = get_stage_vectors('base');
[Q fval] = linprog( c, [], [], A_full, b_vec, LB, UB );

read_results(Q,Period,cellInputFile,solutionFile);

Q = roundn(Q',-2);

%% Check results against original code

disp(' ')
Cost = c';
check_folder = inputLocation;
check_matrix_diffs
