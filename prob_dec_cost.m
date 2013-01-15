function [obj_out] = prob_dec_cost(numscen, dec_cost, gammaprime, mode)

% prob_dec_cost attempts to lower bound the probability of a decrease in
% the optimal value of the LRO.

if nargin < 1
    numscen = [2 4 6 8 10];
end
if nargin < 4
    mode = 'dual';
    if nargin < 3
        gammaprime = 0.42;
        if nargin < 2
            c = [15 7 4 1 2];
            dec_cost = lro_generic( c, numscen );
        end
    end
end

N = sum(numscen);
% The constant inside the log is a number less than 1
Nbar = N*(log(N)-1) - log(gammaprime);
lambda0 = 1;
mu0 = 1;

switch mode
    case 'dual'
        tic
        initGuess = [lambda0; mu0;];
        options = optimset( 'Algorithm','interior-point', 'GradObj','on' );
        % x = fminunc( @(x) opt_obj(x,dec_cost,numscen,N,Nbar), initGuess, ...
        %     options );
        [x obj] = fmincon( @(x) opt_obj(x,dec_cost,numscen,N,Nbar), initGuess, ...
            [], [], [], [], ...
            [0;0], [Inf;Inf], [], options );
        % obj = opt_obj(x,dec_cost,numscen,N,Nbar)
        lambda = x(1);
        mu = x(2);
        pworst = lambda*numscen./(mu+dec_cost);
        % sum(pworst)
        obj_out = -obj;
        toc
    case 'primal'
        tic
        initGuess = numscen/N;
        options = optimset( 'Algorithm','interior-point', ...
            'GradObj','on', 'GradConstr','on' );
        [p obj] = fmincon( @(p) primal_opt_obj(p,dec_cost), initGuess, ...
            [],[],ones(size(initGuess)),1, ...
            zeros(size(initGuess)), ones(size(initGuess)), ...
            @(p) mynonlincon(p,numscen,gammaprime), options );
%         sum(p)
        obj_out = obj;
        toc
    otherwise
        error( 'Choose method: dual or primal' )
end

% i = 1;
% j = 4;
% if dec_cost(i) == dec_cost(j)
%     error('Choose different i,j')
% end
% A = p(i)*numscen(j)/(p(j)*numscen(i));
% mup = (dec_cost(j) - A*dec_cost(i))/(A - 1)
% lambdap = p(i)*(mu + dec_cost(i))/numscen(i)
% opt_obj([lambdap;mup],dec_cost,numscen,N,Nbar)
% 
% l = 0.01:0.01:0.2;
% op = zeros(size(l));
% for ii=1:length(l)
%     op(ii) = -opt_obj([l(ii),mup]',dec_cost,numscen,N,Nbar);
% end
% plot(l,op)

% [lambda, mu]
function [obj deriv] = opt_obj(x,dec_cost,numscen,N,Nbar)
Nl = sum(numscen.*dec_cost);
Nlc = sum(numscen.*(1-dec_cost));
obj = [Nbar 1]*x + N*x(1)*log(x(1)) - ...
    x(1)*(Nl*log(x(2)+1) + Nlc*log(x(2)));
deriv = [Nbar + N*(log(x(1))+1) - Nl*log(x(2)+1) - Nlc*log(x(2)), ...
    1 - x(1)*(Nl/(x(2)+1) + Nlc/x(2)) ];

function [obj deriv] = primal_opt_obj(p,dec_cost)
p = p(:);
obj = dec_cost*p;
deriv = dec_cost;

function [c ceq GC GCeq] = mynonlincon(p,numscen,gammap)
p = p(:);
numscen = numscen(:);
N = sum(numscen);
c = log(gammap) + sum(numscen.*log(numscen/N)) - sum(numscen.*log(p));
ceq = [];
GC = -numscen./p;
GCeq = [];