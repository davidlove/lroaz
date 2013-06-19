function OneStageAZ( inputLocation, numYears, Time_Lag )

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

% The PrintConstraint method of LPModel prints constraints in a human
% readable format.  You can call this function with the name of a user
% whose constraint you'd like to see, as in
disp('-------------------------------------------------------------------')
disp('Printing the SAVSARP constraint by name')
lp.PrintConstraint('SAVSARP')

% With the row number of the constraint you'd like
disp('-------------------------------------------------------------------')
disp('Printing the constraint in row 31')
lp.PrintConstraint(31)

% With a vector or row numbers
disp('-------------------------------------------------------------------')
disp('Printing constraints in 1:3')
lp.PrintConstraint(1:3)

% Or with a cell array of user names
disp('-------------------------------------------------------------------')
disp('Printing constraints for CAVSARP, CAP, and DemP_C')
lp.PrintConstraint({'CAVSARP','CAP','DemP_C'})

% An extra numeric argument allows you to specify a desired time period,
% while an extra argument of the string 'loss' expands out the loss
% variable.  These extra arguments can in either order  The (relative)
% time period of each variable is shown in parentheses.  Variables from 
% previous periods will be omitted when they don't exist.
disp('-------------------------------------------------------------------')
disp('Printing constraints other time periods with expanded loss variables')
lp.PrintConstraint('SAVSARP', 1, 'loss')
lp.PrintConstraint('DemP_C', 'loss', 15)

disp('-------------------------------------------------------------------')