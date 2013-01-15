function gather_lro_data_version2()
%gather_lro_data_version2 loops over lroaz to get lots of data
%   Version 2 begins user lroaz_version15 (or later versions) 
%   that might be useful in gettting rid of the mysterious
%   unboundedness errros.

tStart = tic;

dgp = .5;
% gp = dgp:dgp:1-dgp;
gp = .15:.2:.95;

[exitflag,p1,x1,l1,m1,s1,z1] = lroaz_version15(gp(1));
num = 0;
while exitflag ~= 1 && num < 20
    num = num + 1;
    [exitflag,p1,x1,l1,m1,s1,z1] = lroaz_version15(gp(ii),[x1;l1;m1;0]);
end

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
    disp(' ')
    disp(['gp = ' num2str(gp(ii))])
    [exitflag,p1,x1,l1,m1,s1,z1] = lroaz_version15(gp(ii));
    num = 0;
    while exitflag ~= 1 && num < 20
        num = num + 1;
        [exitflag,p1,x1,l1,m1,s1,z1] = lroaz_version15(gp(ii),[x1;l1;m1;0]);
    end
    pWorst(:,ii) = p1';
    x(:,ii) = x1;
    lambda(:,ii) = l1;
    mu(:,ii) = m1;
    scenCosts(:,ii) =s1;
    zCost(:,ii) = z1;
end

time = toc(tStart);
disp(['Total time: ' num2str(time)])
disp(['Average solution time: ' num2str(time/length(gp))])

save('saved_variables.mat','pWorst','x','lambda','mu','scenCosts','zCost')