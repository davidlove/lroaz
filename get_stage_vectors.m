function [c,A,b,l,u,T] = get_stage_vectors(stage,scenario, ConnectionsFile,cellInputFile,Period,periods1,Time_Lag)

persistent b1 UB1 LB1 Cost1 b2 UB2 LB2 Cost2 A_full_1 A_full_2 tech_matrix

persistent Var_name Nr Nc

persistent Abase A_st A_lag cost

persistent scale

if isempty(b1)
    scale = 1;
    
    [Abase,A_st,A_lag,Var_name,Nr,Nc] = build_A(ConnectionsFile, cellInputFile{1});
    A_full_temp = expand_A(Abase,A_st,A_lag,Period,Time_Lag);
    [r,c] = size(Abase);
    A_full_1 = A_full_temp(1:periods1*r,1:periods1*c);
    A_full_2 = A_full_temp(periods1*r+1:end,periods1*c+1:end);
    tech_matrix = A_full_temp(periods1*r+1:end,1:periods1*c);
    tech_matrix = sparse(tech_matrix);
    [b1,UB1,LB1,Cost1,b2,UB2,LB2,Cost2,cost] = two_stage_input(cellInputFile,Period,periods1);
%     For debugging purposes:
%     remake_A_full = [A_full_1   , zeros(size(A_full_1,1),size(A_full_2,2)); ...
%                      tech_matrix, A_full_2];
%     assert( sum(sum(abs( remake_A_full - A_full_temp ))) == 0 );
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
            error('Output og get_stage_vectors(''names''): Var_Name, Nc, Nr')
        end
        c = Var_name;
        A = Nc;
        b = Nr;
    case 'scale'
        if exist('scenario','var')
            scale = scenario;
        end
        c = scale;
    case 'base'
        if nargout > 4
            error('Output og get_stage_vectors(''base''): Abase, A_st, A_lag, cost')
        end
        c = Abase;
        A = A_st;
        b = A_lag;
        l = cost;
    otherwise
        error('Only stages 1 and 2 possible')
end
