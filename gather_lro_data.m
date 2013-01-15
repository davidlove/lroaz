function gather_lro_data()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

gp = 0.1:0.01:0.95;

[p1,x1,l1,m1,s1,z1] = lroaz_version13(gp(1));

pWorst = zeros(length(p1),length(gp));
x = zeros(length(x1),length(gp));
lambda = zeros(length(l1),length(gp));
mu = zeros(length(m1),length(gp));
scenCosts = zeros(length(s1),length(gp));
zCost = zeros(length(z1),length(gp));

pWorst(:,1) = p1';
x(:,1) = x1;
lambda(:,1) = l1;
mu(:,1) = m1;
scenCosts(:,1) =s1;
zCost(:,1) = z1;

for ii = 2:length(gp)
    disp('')
    disp(['gp = ' num2str(gp(ii))])
    [p1,x1,l1,m1,s1,z1] = lroaz_version14(gp(ii));
    pWorst(:,ii) = p1';
    x(:,ii) = x1;
    lambda(:,ii) = l1;
    mu(:,ii) = m1;
    scenCosts(:,ii) =s1;
    zCost(:,ii) = z1;
end



save('saved_variables.mat','pWorst','x','lambda','mu','scenCosts','zCost')