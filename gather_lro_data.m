function gather_lro_data(saveFileName)
%gather_lro_data_version2 loops over lroaz to get lots of data
%   Version 2 begins user lroaz_version15 (or later versions) 
%   that might be useful in gettting rid of the mysterious
%   unboundedness errros.
%   Version 3 includes the function solve_lroaz, which should 
%   mitigate the possibility of bugs in the solution code.
%   Version 4 uses the two best previous solutions to try to 
%   get around the problems with linprog
%   Version 5 integrates the first run into the loop to 
%   help minimize bugs in this code
%   Version 6 includes counts of cuts and stores time to solution

if nargin < 1 || ~exist(saveFileName,'file')
    newFile = true;
    dgp = .01;
    gp = dgp:dgp:1-dgp;
    if nargin < 1
        saveFileName = 'saved_variables.mat';
    end
    iiSet = 1:length(gp);
else
    newFile = false;
    contLoad = load(saveFileName);
    gp = contLoad.gp;
    pWorst = contLoad.pWorst;
    x = contLoad.x;
    lambda = contLoad.lambda;
    mu = contLoad.mu;
    scenCosts = contLoad.scenCosts;
    zCost = contLoad.zCost;
    likelihood = contLoad.likelihood;
    exitFlags = contLoad.exitFlags;
    numRuns = contLoad.numRuns;
    numCuts = contLoad.numCuts;
    timeRuns = contLoad.timeRuns;

    iiSet = find(exitFlags~=1);
    disp(iiSet)
end
    
tStart = tic;

for ii = iiSet
    disp(' ')
    disp(['gp = ' num2str(gp(ii))])
    timeIndiv = tic;
    [p1,x1,l1,m1,s1,z1,r1,eF,n,c] = solve_lroaz(gp(ii));
    timeIndiv = toc(timeIndiv);
    
    if ii==1 && newFile
        pWorst = zeros(length(p1),length(gp));
        x = zeros(length(x1),length(gp));
        lambda = zeros(length(l1),length(gp));
        mu = zeros(length(m1),length(gp));
        scenCosts = zeros(length(s1),length(gp));
        zCost = zeros(length(z1),length(gp));
        likelihood = zeros(length(r1),length(gp));
        exitFlags = zeros(1,length(gp));
        numRuns = zeros(1,length(gp));
        numCuts = zeros(1,length(gp));
        timeRuns = zeros(1,length(gp));
    end
    
    pWorst(:,ii) = p1';
    x(:,ii) = x1;
    lambda(:,ii) = l1;
    mu(:,ii) = m1;
    scenCosts(:,ii) =s1;
    zCost(:,ii) = z1;
    likelihood(:,ii) = r1;
    exitFlags(ii) = eF;
    numRuns(ii) = n;
    numCuts(ii) = c;
    timeRuns(ii) = timeIndiv;

    save(saveFileName,'gp','pWorst','x','lambda','mu', ...
        'scenCosts','zCost','likelihood', ...
        'exitFlags','numRuns','numCuts','timeRuns')
end

time = toc(tStart);
disp(['Total time: ' num2str(time)])
disp(['Average solution time: ' num2str(time/length(iiSet))])

% save('saved_variables.mat','gp','pWorst','x','lambda','mu', ...
%     'scenCosts','zCost','likelihood','exitFlags','numRuns')


function [p1,x1,l1,m1,s1,z1,r1,exitFlag,num,numCuts] = solve_lroaz(gp)

% zPrev = Inf;
solnBest = [];
soln2Best = [];
num = 1;
totalNum = 10;
exitFlag = -3;
numCuts = 0;

while exitFlag ~= 1 && num <= totalNum
    if num == 1
        clear get_stage_vectors
    end
    disp(' ')
    disp(['gammaprime = ' num2str(gp) ', num = ' num2str(num) ...
        ' of ' num2str(totalNum)])
    disp(' ')
    [exitFlag,p1,solnBestNew,soln2BestNew,s1,z1,r1,cuts1] = lroaz(gp,solnBest);
    if isempty(soln2BestNew)
        solnBest = soln2Best;
    else
        soln2Best = soln2BestNew;
        solnBest = solnBestNew;
    end
%     zPrev = z1;
    numCuts = max(numCuts,cuts1);
    num = num+1;
end
x1 = solnBest(1:end-3);
l1 = solnBest(end-2);
m1 = solnBest(end-1);
num = num-1;
