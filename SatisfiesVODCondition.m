function satisfiesVOD = SatisfiesVODCondition( varargin )
% SatisfiesVODCondition returns a vector indicating which samples satisfy the
% value of data condition
%
% Keyword arguments (*required arguments):
%    obs*: vector of observations numbers of the scenarios
%    p*: vector of worst-case probability distribution
%    phi*: PhiDivergence object or string defining the phi-divergence

% -------------------------------------------------------------------
% Input Parsing
% -------------------------------------------------------------------

requiredArgs = {'obs', 'p', 'phi'};

if mod(length(varargin),2) == 1
    error('Arguments must be key, value pairs')
end

for vv = 1:2:length(varargin)
    key = varargin{vv};
    value = varargin{vv+1};
    switch key
        case 'obs'
            obs = value;
        case 'p'
            p = value;
        case 'phi'
            if isa(value, 'PhiDivergence')
                phi = value;
            elseif ischar(value)
                phi = PhiDivergence( value );
            else
                error('Phi must be PhiDivergence object or string')
            end
        otherwise
            error(['Unknown variable ', key])
    end
end

% -------------------------------------------------------------------
% Default variable assignments
% -------------------------------------------------------------------

% -------------------------------------------------------------------
% Required Variables
% -------------------------------------------------------------------

for aa = requiredArgs
    if ~exist(aa{1}, 'var')
        error([aa{1}, ' is required but not defined'])
    end
end

% -------------------------------------------------------------------
% Value Checking
% -------------------------------------------------------------------

p = p(:);
obs = obs(:);
if size(p) ~= size(obs)
    error('Size of p and obs must match')
end

% -------------------------------------------------------------------
% -------------------------------------------------------------------

N = sum(obs);
q = obs/N;

switch phi.divergence
    case 'burg'
        cost_dec = p./q < N/(N+1);
    case 'chi2'
        const = sum(q.*q./p) + sqrt((N+1)/N);
        cost_dec = const < 2*p./q;
    case 'hellinger'
        const = sum(q.*sqrt(p./q));
        cost_dec = const + sqrt(p./q) < 2*N/(N+1);
    case 'mchi2'
        const = 2 * sum(p.*p./q) - (N+1)^2/N^2;
        cost_dec = (p./q).^2 < const;
    otherwise
        error('No rule for divergence')
end

satisfiesVOD = cost_dec;

end

