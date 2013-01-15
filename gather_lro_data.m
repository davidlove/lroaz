function gather_lro_data_version4()
%gather_lro_data_version2 loops over lroaz to get lots of data
%   Version 2 begins user lroaz_version15 (or later versions) 
%   that might be useful in gettting rid of the mysterious
%   unboundedness errros.
%   Version 3 includes the function solve_lroaz, which should 
%   mitigate the possibility of bugs in the solution code.
%   Version 4 uses the two best previous solutions to try to 
%   get around the problems with linprog

tStart = tic;

dgp = .01;
gp = 0.42:dgp:1-dgp;
% gp = .02;

[p1,x1,l1,m1,s1,z1,r1,eF,n] = solve_lroaz(gp(1));

pWorst = zeros(length(p1),length(gp));
x = zeros(length(x1),length(gp));
lambda = zeros(length(l1),length(gp));
mu = zeros(length(m1),length(gp));
scenCosts = zeros(length(s1),length(gp));
zCost = zeros(length(z1),length(gp));
likelihood = zeros(length(r1),length(gp));
exitFlags = zeros(1,length(gp));
numRuns = zeros(1,length(gp));

pWorst(:,1) = p1';
x(:,1) = x1;
lambda(:,1) = l1;
mu(:,1) = m1;
scenCosts(:,1) =s1;
zCost(:,1) = z1;
likelihood(:,1) = r1;
exitFlags(1) = eF;
numRuns(1) = n;

for ii = 2:length(gp)
    disp(' ')
    disp(['gp = ' num2str(gp(ii))])
    [p1,x1,l1,m1,s1,z1,r1,eF,n] = solve_lroaz(gp(ii));
    pWorst(:,ii) = p1';
    x(:,ii) = x1;
    lambda(:,ii) = l1;
    mu(:,ii) = m1;
    scenCosts(:,ii) =s1;
    zCost(:,ii) = z1;
    likelihood(:,ii) = r1;
    exitFlags(ii) = eF;
    numRuns(ii) = n;

    save('saved_variables.mat','gp','pWorst','x','lambda','mu', ...
        'scenCosts','zCost','likelihood','exitFlags','n')
end

time = toc(tStart);
disp(['Total time: ' num2str(time)])
disp(['Average solution time: ' num2str(time/length(gp))])

save('saved_variables.mat','gp','pWorst','x','lambda','mu', ...
    'scenCosts','zCost','likelihood','exitFlags','numRuns')


function [p1,x1,l1,m1,s1,z1,r1,exitFlag,num] = solve_lroaz(gp)

% zPrev = Inf;
solnBest = [];
soln2Best = [];
num = 1;
totalNum = 10;
exitFlag = -3;

while exitFlag ~= 1 && num <= totalNum
    if num == 1
        clear get_stage_vectors
    end
    disp(' ')
    disp(['gammaprime = ' num2str(gp) ', num = ' num2str(num) ...
        ' of ' num2str(totalNum)])
    disp(' ')
    [exitFlag,p1,solnBestNew,soln2BestNew,s1,z1,r1] = lroaz_version16(gp,solnBest);
    if isempty(soln2BestNew)
        solnBest = soln2Best;
    else
        soln2Best = soln2BestNew;
        solnBest = solnBestNew;
    end
%     zPrev = z1;
    num = num+1;
end
x1 = solnBest(1:end-3);
l1 = solnBest(end-2);
m1 = solnBest(end-1);

% [exitflag,p1,x1,l1,m1,s1,z1,r1] = lroaz_version15(gp);
% zPrev = z1;
% num = 1;
% totalnum = 10;
% while exitflag ~= 1 && num < totalnum
%     num = num + 1;
%     zPrev = z1;
%     disp(' ')
%     disp(['gammaprime = ' num2str(gp) ', num = ' num2str(num) ...
%         ' of ' num2str(totalnum)])
%     disp(' ')
%     [exitflag,p1,x1,l1,m1,s1,z1,r1] = lroaz_version15(gp,[x1;l1;m1;0]);
%     disp(' ')
%     disp(['z nonincreasing: ' num2str(zPrev >= z1) ...
%         ', z decreasing: ' num2str(zPrev > z1) ])
%     disp(' ')
%     if zPrev == z1
%         break
%     end
% end