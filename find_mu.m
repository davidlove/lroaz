function [mu] = find_mu( lambda, numScens, subProbs )

if nargin < 3
    subProbs = ceil(10*rand(30,1));
    if nargin < 2
        numScens = ceil(5*rand(30,1));
        if nargin < 1
            lambda = 1;
        end
    end
end
subProbs = subProbs(:);
numScens = numScens(:);

if length(subProbs) ~= length(numScens)
    error('subProbs and numScens must have same length')
end

[maxProb maxProbIndex] = max(subProbs);
mu = length(subProbs)/2*numScens(maxProbIndex)*lambda + maxProb;
% ymin = mu - sum( lambda * numScens.*log(mu - subProbs) );
% x = (min(subProbs)-1):0.001:(max(subProbs)+length(subProbs)*max(numScens));
% y = zeros(size(x));
% for ii=1:length(y)
%     y(ii) = x(ii) - sum( lambda * numScens.*log(x(ii)-subProbs) );
% end
% plot(x,y,'b', mu,ymin,'ko')
% hold on

% while true
for ii=1:200
%     pause(0.5)
    muOld = mu;
    mu = mu + (sum(lambda * numScens./(mu - subProbs))-1) / ...
        sum(lambda * numScens./((mu - subProbs).^2));
%     ymin = mu - sum( lambda * numScens.*log(mu - subProbs) );
%     plot(x,y,'b', mu,ymin,'ko')
    if abs(mu - muOld) < min(mu,muOld)*0.01
        break
    end
end
mu = max(mu, maxProb+1);
% hold off