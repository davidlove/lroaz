function gather_lro_data_version3()
%gather_lro_data_version2 loops over lroaz to get lots of data
%   Version 2 begins user lroaz_version15 (or later versions) 
%   that might be useful in gettting rid of the mysterious
%   unboundedness errros.
%   Version 3 includes the function solve_lroaz, which should 
%   mitigate the possibility of bugs in the solution code.

tStart = tic;

dgp = .01;
gp = dgp:dgp:1-dgp;
% gp = .15:.2:.95;

[p1,x1,l1,m1,s1,z1,r1] = solve_lroaz(gp(1));

pWorst = zeros(length(p1),length(gp));
x = zeros(length(x1),length(gp));
lambda = zeros(length(l1),length(gp));
mu = zeros(length(m1),length(gp));
scenCosts = zeros(length(s1),length(gp));
zCost = zeros(length(z1),length(gp));
likelihood = zeros(length(r1),length(gp));

pWorst(:,1) = p1';
x(:,1) = x1;
lambda(:,1) = l1;
mu(:,1) = m1;
scenCosts(:,1) =s1;
zCost(:,1) = z1;
likelihood(:,1) = r1;

for ii = 2:length(gp)
    disp(' ')
    disp(['gp = ' num2str(gp(ii))])
    [p1,x1,l1,m1,s1,z1,r1] = solve_lroaz(gp(ii));
    pWorst(:,ii) = p1';
    x(:,ii) = x1;
    lambda(:,ii) = l1;
    mu(:,ii) = m1;
    scenCosts(:,ii) =s1;
    zCost(:,ii) = z1;
    likelihood(:,ii) = r1;

    save('saved_variables.mat','gp','pWorst','x','lambda','mu', ...
        'scenCosts','zCost','likelihood')
end

time = toc(tStart);
disp(['Total time: ' num2str(time)])
disp(['Average solution time: ' num2str(time/length(gp))])

save('saved_variables.mat','gp','pWorst','x','lambda','mu', ...
    'scenCosts','zCost','likelihood')


function [p1,x1,l1,m1,s1,z1,r1] = solve_lroaz(gp)

[exitflag,p1,x1,l1,m1,s1,z1,r1] = lroaz_version15(gp);
% zPrev = z1;
num = 1;
totalnum = 10;
while exitflag ~= 1 && num < totalnum
    num = num + 1;
    zPrev = z1;
    disp(' ')
    disp(['gammaprime = ' num2str(gp) ', num = ' num2str(num) ...
        ' of ' num2str(totalnum)])
    disp(' ')
    [exitflag,p1,x1,l1,m1,s1,z1,r1] = lroaz_version15(gp,[x1;l1;m1;0]);
    disp(' ')
    disp(['z nonincreasing: ' num2str(zPrev >= z1) ...
        ', z decreasing: ' num2str(zPrev > z1) ])
    disp(' ')
    if zPrev == z1
        break
    end
end