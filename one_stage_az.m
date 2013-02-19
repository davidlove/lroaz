function one_stage_az( inputLocation, numYears, Time_Lag )

% one_stage_az demonstrates the correct way to use the LPModel class to
% solve the single-stage problem with time lag

if nargin < 3
    Time_Lag = 1;
    if nargin < 2
        numYears = 41;
        if nargin < 1
            inputLocation = 'SP C/5/';
        end
    end
end

tic;

% The arguments to LPModel are:
%     1. The folder containing the spreadsheets:
%         Connections.xlsx
%         Inputs.xlsx
%         sf.xlsx
%        Here, I've used the folder 'SP C/5'
%     2. The number of years in the model.  This is NOT the number of time
%     periods.
%     3. The time lag.  LPModel figures out timePeriods = numYears *
%     Time_Lag on its own
%     LPModel returns the object lp, which could be given any name
lp = LPModel(inputLocation,numYears,Time_Lag);

% Now, we solve the linear program:
%     min  c*x
%     s.t. A*x = b
%          l <= x <= u
% using the lp object.
% All matrices and vectors are stored within lp
[Q fval] = linprog( lp.c, [], [], lp.A, lp.b, lp.l, lp.u );

% Finally, we read the results given by the solution Q.
% ReadResults also returns the modified (rounded and transposed) version of
% Q
Q = lp.ReadResults(Q,inputLocation);

toc
