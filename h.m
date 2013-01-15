function [fval slope int] = h(x,num)

% h solves the second stage problem of the LRO water implementation.
% Necessary variables:
%   x -- First stage solution
if nargin < 2
    num = 1;
end

[q,D,d,l,u,B] = get_stage_vectors(2,num);

[~,fval,~,~,pi] = linprog(q,[],[],D,d+B*x,l,u);

% Note: I'm pretty sure that Matlab is returning Lagrange multipliers that
% are opposite in sign to those that I calculate.
slope = -pi.eqlin'*B;
int = -pi.eqlin'*d - pi.upper(u<Inf)'*u(u<Inf) - pi.lower(l~=0)'*l(l~=0);
