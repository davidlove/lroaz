function A_full = expand_A(A,Ap,Period)

xc = size(A,1);
xd = size(A,2);
xg = size(Ap,1);
xh = size(Ap,2);

if Period > 1 % condition for solution over 1 year period
    A_full(1:xg,1:xh)=Ap;
    for xe = 2:Period
        xf = xe - 1;
        if xe == Period
            A_full(xf*xc+1:xe*xc,xf*xd+1:xe*xd) = A;
        else
            A_full(xf*xc+1:xf*xc+xg,xf*xd+1:xe*xd) = Ap;
        end
        
    end
else A_full = A;
end