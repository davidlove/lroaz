function Q = one_stage_az( inputLocation )

% clear get_stage_vectors

if nargin < 1
    inputLocation = 'all_scenarios/5/';
end

% one_stage_az uses the code from lroaz_version15 to do a one-stage problem
% that takes a given input

ConnectionsFile = [ inputLocation 'Connections.xlsx'];
cellInputFile = { [inputLocation 'Inputs.xlsx'] };

% Problem Parameters
Period = 41;
periods1 = 41;

[c,A_full,b_vec,LB,UB] = get_stage_vectors(1,1, ...
    ConnectionsFile,cellInputFile,Period,periods1);
[Nc,Nr,Cuser] = get_stage_vectors('names');
[A,Ap,cost] = get_stage_vectors('base');
[Q fval] = linprog( c, [], [], A_full, b_vec, LB, UB );

%% Write Variable Names to Input file

Csource = Nc';
xa = length(Cuser);
xa = Csource(xa+1:end);
Cuser = vertcat(Cuser',xa);
Var_name = cell(length(Cuser),1);
for xb = 1:length(Cuser)
    if ismember(Csource(xb,1),Cuser(xb,1))
        Var_name(xb,1) = Csource(xb,1);
    else
        Var_name(xb,1) = strcat(Csource(xb,1),{' -->> '},Cuser(xb,1));
    end
end

% if Question1 == 1
%     disp('Updating Input file variables...');
%     
%     xlswrite(InputFile,Var_name,'Costs','B2');
%     xlswrite(InputFile,Var_name,'Upper Bounds','B2');
%     xlswrite(InputFile,Var_name,'Lower Bounds','B2');
%     xlswrite(InputFile,Var_name,'kWh','B2');
%     sound(y, Fs);
%     
%     Question2 = input('Enter 1 after "Inputs" file has been updated:  ');
%     clc;
% end

%% User Constraint Check

Q = Q';
Sol_Matrix = zeros(size(A_full));

for xa = 1:size(A_full,1)
    xb = A_full(xa,:);
    xc = xb.*Q;
    Sol_Matrix(xa,:) = xb.*Q;
    xe(xa,1) = sum(xc);
end

Residuals = round(b_vec - xe);
check = any(Residuals);
clc;

%% Sort solution by year

Q = roundn(Q,0);
Flow = zeros(length(A),Period);
Flow(1:length(A),1)= Q(1:length(A));

for xa = 1:Period-1
    xb = length(A)*xa;
    Flow(1:length(A),xa+1) = Q(xb+1:xb+length(A));
end



%% Compute additional results

disp('Computing additional outputs...')

% Individual arc costs (dummy flow set to 0)
blank = cell(250,250);
xa = find(cost(:,1) == 1e10);
cost(xa,:) = 0;
Var_cost = Flow.*cost;

% Total system Losses
xa = strfind(Var_name,'loss');
xa = find(~cellfun(@isempty,xa));
Total_Loss = sum(sum(Flow(xa,:)));

% Total system demand
xa = strfind(Var_name,'Dem');
xa = find(~cellfun(@isempty,xa));
xb = strfind(Var_name,'Dum');
xb = find(~cellfun(@isempty,xb));
xa=setdiff(xa,xb);
xb = strfind(Var_name,'loss');
xb = find(~cellfun(@isempty,xb));
xc=setdiff(xa,xb);
Total_Demand = sum(sum(Flow(xc,:))); % Demand plus end user losses

xa = strfind(Var_name,'Dem');
xa = find(~cellfun(@isempty,xa));
xb = strfind(Var_name,'Dum');
xb = find(~cellfun(@isempty,xb));
xa=setdiff(xa,xb);
xd = setdiff(xa,xc);
Demand_Losses = sum(sum(Flow(xd,:))); % End user losses
Total_Demand = Total_Demand-Demand_Losses;
Total_Demand = Total_Demand + Total_Loss; % Real demand is initial demand plus system losses
Percent_Loss = Total_Loss/Total_Demand;

% Total Demand Shortages
xa = strfind(Var_name,'Dum');
xa = find(~cellfun(@isempty,xa));
Total_Shortages = sum(sum(Flow(xa,:)));

% Total Potable Storage
xa = strfind(Var_name,'strg');
xa = find(~cellfun(@isempty,xa));
xb = strfind(Var_name,'NPstrg');
xb = find(~cellfun(@isempty,xb));
xa = setdiff(xa,xb);
Potable_Storage = sum(sum(Flow(xa,end)));

% Energy Consumption
% Energy = zeros(size(Flow));
% for xa = 1:length(kWh)
%     for xb = 1:size(Flow,2)
%         Energy(xa,xb) = kWh(xa)*Flow(xa,xb);
%     end
% end
% Total_Energy = sum(sum(Energy));

% GHG Emissions
% GHG = 0.000623854*Total_Energy;

% Total Reuse
xa = strfind(Var_name,'DemNP');
xa = find(~cellfun(@isempty,xa));
xb = strfind(Var_name,'Dum');
xb = find(~cellfun(@isempty,xb));
xa=setdiff(xa,xb);
xb = strfind(Var_name,'loss');
xb = find(~cellfun(@isempty,xb));
xc=setdiff(xa,xb);
xa = strmatch('P',Var_name);
xb = setdiff(xc,xa);
Total_Reuse = sum(sum(Flow(xb,:)));
Percent_Reuse = Total_Reuse/Total_Demand;


%% Write to excel solution file

disp('Writing to excel...')

Flow = num2cell(Flow);
Flow = horzcat(Var_name,Flow);
cost = num2cell(cost);
cost = horzcat(Var_name,cost);
%     cost = cost(:,1:Period);
Total_cost = sum(sum(Var_cost));
% Final = [Total_Demand,Total_Loss,Percent_Loss,Total_Shortages,...
%     Potable_Storage,Total_cost/1e6,...
%     GHG/1e6,Total_Reuse,Percent_Reuse];
Var_cost = num2cell(Var_cost);
Var_cost = horzcat(Var_name,Var_cost);
% if Question1 == 1
%     xlswrite('Solution.xlsx',blank,'Solution','A2');
%     xlswrite('Solution.xlsx',Flow,'Solution','A2');
%     xlswrite('Solution.xlsx',blank,'Costs','A2');
%     xlswrite('Solution.xlsx',cost,'Costs','A2');
%     xlswrite('Solution.xlsx',blank,'Price','A2');
%     xlswrite('Solution.xlsx',Var_cost,'Price','A2');
%     xlswrite('Solution.xlsx',Final,'Final Stats','B5');
% else
%     xlswrite('Solution.xlsx',Flow,'Solution','A2');
%     xlswrite('Solution.xlsx',cost,'Costs','A2');
%     xlswrite('Solution.xlsx',Var_cost,'Price','A2');
%     xlswrite('Solution.xlsx',Final,'Final Stats','B5');
% end

% toc;
% xf = toc;
clc;


%% System details

disp(' ');
disp(' =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = ');
disp('                                      Sytem Outline');
disp(' =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  =  = ');

% disp(['Time to load data/construct matrix(min)  =  ',num2str(xf/60)]);
disp(' ');
% disp(['Number of zones  =  ',num2str(floor(max(Zones)))]);
% disp(' ');
disp(['Number of arcs  =  ',num2str(length(A))]);
disp(' ');
disp(['Period (yr)  =  ',num2str(Period)]);

disp(' ')
Cost = c';
check_folder = inputLocation;
check_matrix_diffs
