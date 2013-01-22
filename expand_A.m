function A_full = expand_A(A,A_st,A_lag,Period,Time_Lag)

xc = size(A,1); % Initial Matrix without storage carry over
xd = size(A,2);
% xg = size(A_st,1); % Initial Matrix with storage carry over (incorporates two time periods)
% xh = size(A_st,2);
A_full = zeros(xc*Period,xd*Period);
A_full(1:xc,1:xd)=A;
A_lag_ini = A_st + A_lag;
if Period > 1 % condition for more than one time step
    for xe = 2:Period
        xf = xe - 1;
        A_full(xf*xc+1:xe*xc,xf*xd+1:xe*xd)=A;
        if Period/Time_Lag == Period
            A_full(xf*xc+1:xe*xc,(xf-1)*xd+1:xf*xd)=A_lag_ini;
        elseif xf >= Time_Lag
            A_full(xf*xc+1:xe*xc,(xf-Time_Lag)*xd+1:(xf-Time_Lag+1)*xd)=A_lag;
            A_full(xf*xc+1:xe*xc,(xf-1)*xd+1:xf*xd)=A_st;
        else
            A_full(xf*xc+1:xe*xc,(xf-1)*xd+1:xf*xd)=A_st;
        end
    end
end