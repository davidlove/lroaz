function [wwtpUsage, iprUsage] = MaxInfrastructureUsage( savedDataFile, lp )
%MaxInfrastructureUsage computes the necessary capacity, in MGD, for both
%the satellite WWTP and decentralized IPR

d = load(savedDataFile);
[scens, rhos] = size(d.scenSolns);
numVars = size(lp.Abase, 2);
secondStagePeriods = length(d.scenSolns{1,1}) / numVars;
wwtpLocs = find(~cellfun(@isempty, regexp(lp.variableNames, '-->> SP_C')));
iprLocs = find(~cellfun(@isempty, regexp(lp.variableNames, '-->> IPR_C')));

wwtpUsage = zeros(scens, rhos);
iprUsage = zeros(scens, rhos);

for ii=1:scens
    for jj=1:rhos
        solns = d.scenSolns{ii,jj};
        solns = reshape(solns, numVars, secondStagePeriods);
        wwtpUsage(ii,jj) = max(sum(solns(wwtpLocs,:), 1));
        iprUsage(ii,jj) = max(sum(solns(iprLocs,:), 1));
    end
end

wwtpUsage = max(max(wwtpUsage)) / 1121; % Convert to MGD
iprUsage = max(max(iprUsage)) / 1121; % Convert to MGD

end
