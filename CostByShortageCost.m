function newCost = CostByShortageCost( savedDataFile, lp, shortageCost, rho )

%CostByShortageCost changes the shortage cost and recalculates the expected
%worst-case cost for any solution

d = load(savedDataFile);
[scens, rhos] = size(d.scenSolns);
numVars = size(lp.Abase, 2);
firstStagePeriods = length(lp.c) / numVars;
secondStagePeriods = length(d.scenSolns{1,1}) / numVars;

dummyLocs = find(~cellfun(@isempty, regexp(lp.variableNames, 'Dummy -->>')));
c = lp.c;
for ii=1:firstStagePeriods
    locs = dummyLocs + (ii-1)*numVars;
    c(locs) = shortageCost / ((1+0.04)^(ii-1));
end
firstStageCost = c*d.x(:,rho);

secondStageCosts = zeros(scens, 1);
for omega=1:scens
    q = lp.Getq(omega);
    for ii = 1:secondStagePeriods
        locs = dummyLocs + (ii-1)*numVars;
        q(locs) = shortageCost / ((1+0.04)^(firstStagePeriods+ii-1));
    end
    secondStageCosts(omega) = q*d.scenSolns{omega,rho};
end

pWorst = d.pWorst(:,rho)';
pWorst = pWorst / sum(pWorst);
newCost = firstStageCost + pWorst * secondStageCosts;
