function A_full = expand_A(A,A_st,A_lag,Period,Time_Lag)

[rowsA,colsA] = size(A);

if Time_Lag == 1
    A_st = A_st + A_lag;
    assert(issparse(A_st));
else
    [iA_lag_full,jA_lag_full,sA_lag_full] = GenerateDiagonal( A_lag, Time_Lag, Period);
end

[iA_full,jA_full,sA_full] = GenerateDiagonal( A, 0, Period );
[iA_st_full,jA_st_full,sA_st_full] = GenerateDiagonal( A_st, 1, Period);

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
% assert(sum(sum(A_full~=A_full_old))==0);

function [iOut,jOut,sOut] = GenerateDiagonal( matrixIn, subDiag, Period )

% GenerateDiagonal generates the sparse matrix coordinates of the block
% diagonals and sub-diagonals of the constraint matrix A_full. It generates
% copies of matrixIn on the subDiag'th sub-diagonal of a matrix with Period
% periods.
% 
% Note: subDiag = 0 is the main block diagonal.
%       subDiag = 1 is the first block sub-diagonal, etc.

assert(subDiag >= 0)
assert(subDiag < Period)

[rowsM,colsM] = size(matrixIn);

[iIn, jIn, sIn]  = find(matrixIn);
nzM = length(sIn);

iOut = zeros( (Period-subDiag)*length(sIn), 1 );
jOut = iOut;
sOut = iOut;

for pp=subDiag+1:Period
    iOut( (pp-subDiag-1)*nzM + (1:nzM) ) = iIn + (pp-1)*rowsM;
    jOut( (pp-subDiag-1)*nzM + (1:nzM) ) = jIn + (pp-subDiag-1)*colsM;
    sOut( (pp-subDiag-1)*nzM + (1:nzM) ) = sIn;
end