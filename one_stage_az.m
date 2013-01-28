function Q = one_stage_az( inputLocation, Period, Time_Lag )

clear get_stage_vectors
clc

if nargin < 3
    Time_Lag = 1;
    if nargin < 2
        Period = 41;
        if nargin < 1
            inputLocation = 'all_scenarios_new/9/';
        end
    end
end

% one_stage_az uses the code from lroaz_version15 to do a one-stage problem
% that takes a given input

tic;


lp = LPModel(inputLocation,Period,Time_Lag);
[Q fval] = linprog( lp.c, [], [], lp.A, lp.b, lp.l, lp.u );
lp.ReadResults(Q,inputLocation);

toc