function read_results(Q,cost,Period)

%% Get needed informationn

[Nc,~,Cuser] = get_stage_vectors(-1);
Flow = reshape(roundn(Q,0),175,Period);
cost = reshape(cost,size(Flow));

%% Alicia's code to calculate Var_name

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

%% Alicia's code to find various things of interest 
 
 % Individual arc costs (dummy flow set to 0)
% blank = cell(250,250);
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

% % Energy Consumption
% Energy = zeros(size(Flow));
% for xa = 1:length(kWh)
%     for xb = 1:size(Flow,2)
%         Energy(xa,xb) = kWh(xa)*Flow(xa,xb);
%     end
% end
% Total_Energy = sum(sum(Energy));
% 
% % GHG Emissions
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

disp(['Percent Loss = ' num2str(Percent_Loss)])
disp(['Total Shortage = ' num2str(Total_Shortages)])
disp(['Potable Storage = ' num2str(Potable_Storage)])
disp(['Percent Reuse = ' num2str(Percent_Reuse)])