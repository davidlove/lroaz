% function [ demCoeffs ] = WaterDemandModel( spreadsheet )
%WaterDemandModel builds a prediction for water demand from historical data

[pop,~] = xlsread(spreadsheet, 'Population');
[services,servtxt] = xlsread(spreadsheet, 'Services');
[demand,demtxt] = xlsread(spreadsheet, 'Sheet1');

popPredictors = 3:6;
disp(' ')
disp('Building population prediction model based on predictors:')
disp(strjoin(servtxt(popPredictors), ', '))
disp(' ')

basemonth = 7;
% for basemonth = 1:12
popCoeffs = services(basemonth:12:end,popPredictors) \ pop(:,2);
% errors = mean(abs(pop(:,2) - services(basemonth:12:end,popPredictors)*popCoeffs));
% disp([basemonth, errors])
% end

demPredictors = 4:5;
demand = [demand, services(:,popPredictors)*popCoeffs];
demtxt = [demtxt, {'Population'}];

disp('Building water demand model based on predictors:')
disp(strjoin(demtxt(demPredictors), ', '))
disp(' ')

demLocation = 3;
assert(isequal( demtxt{demLocation}, 'Demand' ))

daysPerMonth = [31,28,31,30,31,30,31,31,30,31,30,31]';
daysPerMonth = repmat(daysPerMonth, length(pop(:,2)), 1);
demand(:,demLocation) = demand(:,demLocation) ./ daysPerMonth * 1e6;

demCoeffs = [ones(size(demand,1),1), demand(:,demPredictors)] ...
    \ (demand(:,demLocation)./demand(:,end));
