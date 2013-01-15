clear all
clc

% DebugStageConnection attempts to figure out why Matlab says that some of
% the subproblems are infeasible.

ConnectionsFile = 'all_scenarios/5/Connections.xlsx';
cellInputFile = { ...
    'all_scenarios/5/Inputs.xlsx', ...
    'all_scenarios/6/Inputs.xlsx', ...
    'all_scenarios/7/Inputs.xlsx', ...
    'all_scenarios/8/Inputs.xlsx', ...
    };

% Problem Parameters
numscen = 10*ones(size(cellInputFile));
N = sum(numscen);
% The constant inside the log is a number less than 1
gammaprime = 0.5;
Nbar = N*(log(N)-1) - log(gammaprime);
Period = 4;
periods1 = 2;
% c = 1;

[c,A,Rhs,l,u] = get_stage_vectors(1,1, ...
    ConnectionsFile,cellInputFile,Period,periods1);
x0 = linprog( c, [], [], A, Rhs, l, u );
A = [A zeros(size(A,1),3)];

% Initialize everything
% x0 = -1;
lambda0 = 1;
zlower = -Inf;
zupper = Inf;
objA = [];
objRhs = [];
feasA = [];
feasRhs = [];
feasSlope = [];
feasInt = [];

% Linear options
% options = optimset('MaxIter',85);

% Nonlinear options
options = optimset('MaxIter',1000, ...
    'Algorithm','interior-point', ...
    'GradObj','on');

% Uniform lower bounds on scenario costs
scenLowBnd = 0;



% ------------------------------------------------------------
% -------------  SCRIPT STUFF --------------------------------
% ------------------------------------------------------------

[Nc,Nr,Cuser] = get_stage_vectors(-1);
[q1,D1,d1,l1,u1,B1] = get_stage_vectors(2,1);
[q3,D3,d3,l3,u3,B3] = get_stage_vectors(2,3);

for ii=1:length(numscen)
    h(x0,ii);
end