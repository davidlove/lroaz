function [b1,UB1,LB1,Cost1,b2,UB2,LB2,Cost2,costOut] = BuildVectors(obj, cellInputFile)

% BuildVectors builds all vectors (right hand side b, upper bound UB, lower
% bound LB, cost vector Cost) for both the first and second stages.  First
% stage vectors names have a 1 at the end, second stage vector names have a
% 2 at the end.

numFiles = length(cellInputFile);

for ii=1:numFiles
    
    % Unit costs
    c = xlsread(cellInputFile{ii},'Costs');
    costVector(:,:,ii) = c(2:end,:);
    
    % Variable Bounds
    u = xlsread(cellInputFile{ii},'Upper Bounds');
    ub(:,:,ii) = u(2:end,:);
    l = xlsread(cellInputFile{ii},'Lower Bounds');
    lb(:,:,ii) = l(2:end,:);
    
    % b vector (knowns)
    R = xlsread(cellInputFile{ii},'b_vec');
    RHS(:,:,ii) = R(2:end,:);
    
end

f = find(ub==-1);
ub(f) = Inf(size(f));

if obj.timeLag ~= 1
    assert(obj.timeLag == 12);
    
    numVariablesPerPeriod = size(costVector,1);
    numConstraintsPerPeriod = size(RHS,1);
    
    assert( size(lb,1) == numVariablesPerPeriod );
    assert( size(ub,1) == numVariablesPerPeriod );
    
    costExp = zeros( numVariablesPerPeriod, obj.timePeriods, numFiles );
    rhsExp = zeros( numConstraintsPerPeriod, obj.timePeriods, numFiles );
    lbExp = zeros( numVariablesPerPeriod, obj.timePeriods, numFiles );
    ubExp = zeros( numVariablesPerPeriod, obj.timePeriods, numFiles );
    
    for ff = 1:numFiles
        sf = xlsread( strrep(cellInputFile{ff},'Inputs','sf'),'sf' );
        
        for jj = 1:obj.numYears
            for ii = 1:obj.timeLag
                kk = (jj-1)*obj.timeLag + ii;
                costExp(:,kk,ff) = costVector(:,jj,ff);
                ubExp(:,kk,ff) = ub(:,jj,ff) / obj.timeLag;
                lbExp(:,kk,ff) = lb(:,jj,ff) / obj.timeLag;
                
                rhsExp(1:15,kk,ff) = RHS(1:15,jj,ff) / obj.timeLag;
                rhsExp(16:end,kk,ff) = RHS(16:end,jj,ff) .* sf(ii,1);
                rhsExp(64:100,kk,ff) = RHS(64:100,jj,ff) / obj.timeLag;
            end
        end
        costVector = costExp;
        RHS = rhsExp;
        lb = lbExp;
        ub = ubExp;
    end
end

% if max(abs([size(cost,2),size(ub,2),size(lb,2),size(RHS,2)]-obj.timePeriods)) > 0
if obj.timePeriods > min([size(costVector,2),size(ub,2),size(lb,2),size(RHS,2)])
    error('Request for more periods than data provides');
end

costVector = costVector(:,1:obj.timePeriods,:);
ub = ub(:,1:obj.timePeriods,:);
lb = lb(:,1:obj.timePeriods,:);
RHS = RHS(:,1:obj.timePeriods,:);

% B_vec = RHS; % Make copy for post analysis purposes

% Augment the initial storage by the required amounts
% for ii=1:length(st_nodes)
%      RHS(strcmp(Nr,st_nodes{ii}),1) = RHS(strcmp(Nr,st_nodes{ii}),1) ...
%          - st_initial(ii);
% end

for ii=1:numFiles
    b1t = RHS(:,1:obj.firstStagePeriods,ii);
    b1(:,ii) = b1t(:);
    b2t = RHS(:,obj.firstStagePeriods+1:end,ii);
    if ii > 1
        for jj=1:size(b2t,2)
            b2t(:,jj) = b2t(:,jj) ...
                - RHS(:,obj.firstStagePeriods,ii) + RHS(:,obj.firstStagePeriods,1);
        end
    end
    b2(:,ii) = b2t(:);
    
    
    UB1t = ub(:,1:obj.firstStagePeriods,ii);
    UB1(:,ii) = UB1t(:);
    UB2t = ub(:,obj.firstStagePeriods+1:end,ii);
    if ii > 1
        index = ub(:,1,1) < Inf;
        for jj=1:size(UB2t,2)
            UB2t(index,jj) = UB2t(index,jj) ...
                - ub(index,obj.firstStagePeriods,ii) + ub(index,obj.firstStagePeriods,1);
        end
    end
    UB2(:,ii) = UB2t(:);
    
    LB1t = lb(:,1:obj.firstStagePeriods,ii);
    LB1(:,ii) = LB1t(:);
    LB2t = lb(:,obj.firstStagePeriods+1:end,ii);
    LB2(:,ii) = LB2t(:);
    
    Cost1t = costVector(:,1:obj.firstStagePeriods,ii);
    Cost1(:,ii) = Cost1t(:);
    Cost2t = costVector(:,obj.firstStagePeriods+1:end,ii);
    Cost2(:,ii) = Cost2t(:);
end

b1 = b1(:,1);
UB1 = UB1(:,1);
LB1 = LB1(:,1);
Cost1 = Cost1(:,1);

% b_vec = RHS(:);
% UB = ub(:);
% LB = lb(:);
% Cost = cost(:);

costOut = costVector(:,:,1);
