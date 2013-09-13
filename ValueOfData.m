function [predDecrease, actualDecrease] = ValueOfData( lp, phi, obs, rho )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N = sum(obs);
q = obs(:)/N;

l = length(rho);
n = length(obs);

predDecrease = zeros(n,l);
actualDecrease = predDecrease;

for ii=1:l
    plpBase = SolveLRLP('lp',lp, 'phi',phi, 'obs',obs, 'rho',rho(ii));
    p = plpBase.pWorst;
    p = p(:);
    
    switch phi.divergence
        case 'burg'
            cost_dec = p./q < N/(N+1);
        case 'mchi2'
            const = 2 * sum(p.*p./q) - (N+1)^2/N^2;
            cost_dec = (p./q).^2 < const;
        otherwise
            error('No rule for divergence')
    end
    predDecrease(:,ii) = cost_dec(:);
    for jj=1:n
        newObs = obs;
        newObs(jj) = newObs(jj) + 1;
        plpMod = SolveLRLP('lp',lp, 'phi',phi, ...
            'obs',newObs, 'rho',N/(N+1)*rho(ii));
        actualDecrease(jj,ii) = 1 - plpMod.ObjectiveValue/plpBase.ObjectiveValue;
        if predDecrease(jj,ii) && ~actualDecrease(jj,ii)
            disp(['No cost decrease, (i,j) = (' num2str(ii), ',' num2str(jj), ')'])
        end
    end
end

end

