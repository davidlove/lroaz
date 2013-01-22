function A_full = expand_A(A,A_st,A_lag,Period,Time_Lag)

rowsA = size(A,1); % Initial Matrix without storage carry over
colsA = size(A,2);
% xg = size(A_st,1); % Initial Matrix with storage carry over (incorporates two time periods)
% xh = size(A_st,2);
A_full = zeros(rowsA*Period,colsA*Period);
A_full(1:rowsA,1:colsA)=A;
A_lag_ini = A_st + A_lag;
% if Period > 1 % condition for more than one time step
for xe = 2:Period
    xf = xe - 1;
    A_full(xf*rowsA+1:xe*rowsA,xf*colsA+1:xe*colsA)=A;
%     if Period/Time_Lag == Period
    if Time_Lag == 1
        A_full(xf*rowsA+1:xe*rowsA,(xf-1)*colsA+1:xf*colsA)=A_lag_ini;
    elseif xf >= Time_Lag
        A_full(xf*rowsA+1:xe*rowsA,(xf-Time_Lag)*colsA+1:(xf-Time_Lag+1)*colsA)=A_lag;
        A_full(xf*rowsA+1:xe*rowsA,(xf-1)*colsA+1:xf*colsA)=A_st;
    else
        A_full(xf*rowsA+1:xe*rowsA,(xf-1)*colsA+1:xf*colsA)=A_st;
    end
end
% end