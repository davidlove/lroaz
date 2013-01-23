function A_full = expand_A(A,A_st,A_lag,Period,Time_Lag)

[rowsA,colsA] = size(A);
% rowsA = size(A,1); % Initial Matrix without storage carry over
% colsA = size(A,2);
% xg = size(A_st,1); % Initial Matrix with storage carry over (incorporates two time periods)
% xh = size(A_st,2);
A_full_old = zeros(rowsA*Period,colsA*Period);

A_full_old(1:rowsA,1:colsA)=A;
A_lag_ini = A_st + A_lag;
% if Period > 1 % condition for more than one time step
for xe = 2:Period
    xf = xe - 1;
    A_full_old(xf*rowsA+1:xe*rowsA,xf*colsA+1:xe*colsA)=A;
    %     if Period/Time_Lag == Period
    if Time_Lag == 1
        A_full_old(xf*rowsA+1:xe*rowsA,(xf-1)*colsA+1:xf*colsA)=A_lag_ini;
    elseif xf >= Time_Lag
        A_full_old(xf*rowsA+1:xe*rowsA,(xf-Time_Lag)*colsA+1:(xf-Time_Lag+1)*colsA)=A_lag;
        A_full_old(xf*rowsA+1:xe*rowsA,(xf-1)*colsA+1:xf*colsA)=A_st;
    else
        A_full_old(xf*rowsA+1:xe*rowsA,(xf-1)*colsA+1:xf*colsA)=A_st;
    end
end
% end

if Time_Lag == 1
    A_st = A_st + A_lag;
    assert(issparse(A_st));
else
    % Do nothing
end

[iA,    jA,    sA]     = find(A);
[iA_st, jA_st, sA_st]  = find(A_st);
[iA_lag,jA_lag,sA_lag] = find(A_lag);

iA_full = zeros(Period*length(sA),1);
jA_full = iA_full;
sA_full = iA_full;

iA_st_full = zeros( (Period-1)*length(sA_st), 1 );
jA_st_full = iA_st_full;
sA_st_full = iA_st_full;

if Time_Lag == 1
    % Do nothing
else
    iA_lag_full = zeros( (Period-Time_Lag)*length(sA_log), 1 );
    jA_lag_full = iA_lag_full;
    sA_lag_full = iA_lag_full;
end

nzA     = length(sA);
nzA_st  = length(sA_st);
if Time_Lag == 1
    % Do nothing
else
    nzA_lag = length(sA_lag);
end

for pp=1:Period
    iA_full( (pp-1)*nzA + (1:nzA) ) = iA + (pp-1)*rowsA;
    jA_full( (pp-1)*nzA + (1:nzA) ) = jA + (pp-1)*colsA;
    sA_full( (pp-1)*nzA + (1:nzA) ) = sA;
end

for pp=2:Period
    iA_st_full( (pp-2)*nzA_st + (1:nzA_st) ) = iA_st + (pp-1)*rowsA;
    jA_st_full( (pp-2)*nzA_st + (1:nzA_st) ) = jA_st + (pp-2)*colsA;
    sA_st_full( (pp-2)*nzA_st + (1:nzA_st) ) = sA_st;
end

if Time_Lag == 1
    % Do nothing
else
    for pp=(Time_Lag+1):Period
        
    end
end

if Time_Lag == 1
    i = [iA_full; iA_st_full];
    j = [jA_full; jA_st_full];
    s = [sA_full; sA_st_full];
else
    i = [iA_full; iA_st_full; iA_lag_full];
    j = [jA_full; jA_st_full; jA_lag_full];
    s = [sA_full; sA_st_full; sA_lag_full];
end

A_full = sparse(i,j,s,Period*rowsA,Period*colsA,length(i));
assert(sum(sum(A_full~=A_full_old))==0);
