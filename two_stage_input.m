function [b1,UB1,LB1,Cost1,b2,UB2,LB2,Cost2,costOut] = two_stage_input(cellInputFile,Period,periods1)

% Reads data from InputFile, splits into vectors for the first and second
% stages of the problem.

if ~iscell(cellInputFile)
    error('Input files must be specified in cell arrays')
end

numcell = length(cellInputFile);

for ii=1:numcell
    
    % Unit costs
    c = xlsread(cellInputFile{ii},'Costs');
    cost(:,:,ii) = c(2:end,:);
    
    % Variable Bounds
    u = xlsread(cellInputFile{ii},'Upper Bounds');
    ub(:,:,ii) = u(2:end,:);
    l = xlsread(cellInputFile{ii},'Lower Bounds');
    lb(:,:,ii) = l(2:end,:);
    
    % b vector (knowns)
    R = xlsread(cellInputFile{ii},'b_vec');
    RHS(:,:,ii) = R(2:end,:);
    
end

% if max(abs([size(cost,2),size(ub,2),size(lb,2),size(RHS,2)]-Period)) > 0
if Period > min([size(cost,2),size(ub,2),size(lb,2),size(RHS,2)])
    error('Request for more periods than data provides');
end

f = find(ub==-1);
ub(f) = Inf(size(f));

cost = cost(:,1:Period,:);
ub = ub(:,1:Period,:);
lb = lb(:,1:Period,:);
RHS = RHS(:,1:Period,:);

% B_vec = RHS; % Make copy for post analysis purposes

% Augment the initial storage by the required amounts
% for ii=1:length(st_nodes)
%      RHS(strcmp(Nr,st_nodes{ii}),1) = RHS(strcmp(Nr,st_nodes{ii}),1) ...
%          - st_initial(ii);
% end

for ii=1:numcell
    b1t = RHS(:,1:periods1,ii);
    b1(:,ii) = b1t(:);
    b2t = RHS(:,periods1+1:end,ii);
    if ii > 1
        for jj=1:size(b2t,2)
            b2t(:,jj) = b2t(:,jj) ...
                - RHS(:,periods1,ii) + RHS(:,periods1,1);
        end
    end
    b2(:,ii) = b2t(:);
    

    UB1t = ub(:,1:periods1,ii);
    UB1(:,ii) = UB1t(:);
    UB2t = ub(:,periods1+1:end,ii);
    if ii > 1
        index = ub(:,1,1) < Inf;
        for jj=1:size(UB2t,2)
            UB2t(index,jj) = UB2t(index,jj) ...
                - ub(index,periods1,ii) + ub(index,periods1,1);
        end
    end
    UB2(:,ii) = UB2t(:);
    
    LB1t = lb(:,1:periods1,ii);
    LB1(:,ii) = LB1t(:);
    LB2t = lb(:,periods1+1:end,ii);
    LB2(:,ii) = LB2t(:);
    
    Cost1t = cost(:,1:periods1,ii);
    Cost1(:,ii) = Cost1t(:);
    Cost2t = cost(:,periods1+1:end,ii);
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

costOut = cost(:,:,1);
