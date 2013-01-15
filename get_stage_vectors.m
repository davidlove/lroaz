function [c,A,b,l,u,T] = get_stage_vectors(stage,scenario, ConnectionsFile,cellInputFile,Period,periods1)

persistent b1 UB1 LB1 Cost1 b2 UB2 LB2 Cost2 A_full_1 A_full_2 tech_matrix

persistent Nr Nc Cuser

persistent Abase Ap cost

persistent scale

if isempty(b1)
    scale = 1;
    
    [Abase,Ap,Nr,Nc,Zones,Cuser] = build_A(ConnectionsFile, cellInputFile{1});
    A_full_1 = expand_A(Abase,Ap,periods1);
    A_full_2 = expand_A(Abase,Ap,Period-periods1);
    cols_strg = find(~cellfun(@isempty,strfind(Nc,'-strg')));
    len_strg = cell2mat(strfind(Nc,'-strg'))-1;
    rows_strg = zeros(size(len_strg));
    for ii=1:length(rows_strg)
         rows_strg(ii) = find(~cellfun(@isempty,strfind(Nr,Nc{cols_strg(ii)}(1:len_strg(ii)))),1,'first');
    end
    tech_matrix = zeros(size(A_full_2,1),size(A_full_1,2));
    tech_matrix(1:size(Abase,1),(periods1-1)*size(Abase,2)+1:end) = ...
        - Ap(size(Abase,1)+1:end,:);
    tech_matrix = sparse(tech_matrix);
    [b1,UB1,LB1,Cost1,b2,UB2,LB2,Cost2,cost] = two_stage_input(cellInputFile,Period,periods1);
end

switch stage
    case 1
        b = b1;
        u = UB1;
        l = LB1;
        c = scale*Cost1';
        A = A_full_1;
    case 2
        if scenario > size(b2,2)
            error('Not enough scenarios')
        end
        b = b2(:,scenario);
        u = UB2(:,scenario);
        l = LB2(:,scenario);
        c = scale*Cost2(:,scenario)';
        A = A_full_2;
        T = tech_matrix;
    case {-1, 'names'}
        if nargout > 3
            error('Can only output Nc, Nr and Cuser in this mode')
        end
        c = Nc;
        A = Nr;
        b = Cuser;
    case 'scale'
        if exist('scenario','var')
            scale = scenario;
        end
        c = scale;
    case 'base'
        c = Abase;
        A = Ap;
        b = cost;
    otherwise
        error('Only stages 1 and 2 possible')
end
