classdef PhiDivergence < handle
    %PhiDivergence defines a phi-divergence for use in robust optimization
    %   Detailed explanation goes here
    
    properties (GetAccess=public, SetAccess=immutable)
        divergence
        limit
    end
    
    properties (GetAccess=private, SetAccess=immutable)
        func
        conjugate
        conjugateDerivative
        phi2Derivative
    end
    
    properties (GetAccess=public, SetAccess=private)
        computationLimit
    end
    
    methods
        function obj = PhiDivergence( inDivergence )
            % Constructor defines which phi-divergence is to be used.
            switch lower(inDivergence)
                case 'lro'
                    obj.func = @(t) -log(t);
                    obj.conjugate = @(s) -log(-s)-1;
                    obj.conjugateDerivative = @(s) -1./s;
                    obj.phi2Derivative = @(t) 1./(t.^2);
                    obj.limit = 0;
                case 'burg'
                    obj.func = @(t) -log(t) + t - 1;
                    obj.conjugate = @(s) -log(1-s);
                    obj.conjugateDerivative = @(s) 1./(1-s);
                    obj.phi2Derivative = @(t) 1./(t.^2);
                    obj.limit = 1;
                case 'kl'
                    % include realmin to prevent 0*log(0) from evaluating
                    % as NaN.
                    obj.func = @(t) t.*log(max(t,realmin)) - t + 1;
                    obj.conjugate = @(s) exp(s) - 1;
                    obj.conjugateDerivative = @(s) exp(s);
                    obj.phi2Derivative = @(t) 1./t;
                    obj.limit = Inf;
                case 'chi2'
                    % abs(t) prevents getting -Inf when t == 0
                    obj.func = @(t) (t-1).^2./abs(t);
                    obj.conjugate = @(s) 2 - 2*sqrt(1-s);
                    obj.conjugateDerivative = @(s) 1./sqrt(1-s);
                    obj.phi2Derivative = @(t) 2./t^3;
                    obj.limit = 1;
                case 'mchi2'
                    obj.func = @(t) (t-1).^2;
                    obj.conjugate = @(s) (-1).*(s < -2) + ...
                        (max(s,-2)+max(s,-2).^2/4).*(s >= -2);
                    obj.conjugateDerivative = @(s) (1+max(s,-2)/2).*(s >= -2);
                    obj.phi2Derivative = @(t) 2;
                    obj.limit = Inf;
                case 'hellinger'
                    obj.func = @(t) (sqrt(t)-1).^2;
                    obj.conjugate = @(s) max(s,-realmax)./(1-max(s,-realmax));
                    obj.conjugateDerivative = @(s) 1./((s-1).^2);
                    obj.phi2Derivative = @(t) 1./(2*t^(3/2));
                    obj.limit = 1;
                otherwise
                    error('PhiDivergence:PhiDivergence:Unknown', ...
                        ['Unknown phi type ' inDivergence])
            end
            obj.divergence = inDivergence;
            obj.computationLimit = Inf;
        end
        
        function outVal = Contribution( obj, inNumer, inDenom )
            zD = inDenom == 0;
            zN = zD & (inNumer == 0);
            outVal = inDenom.*obj.func(inNumer./inDenom);
            outVal(zD) = inNumer(zD) * obj.limit;
            outVal(zN) = 0;
        end
        
        function outVal = Conjugate( obj, inS )
            % Conjugate returns the value of the conjugate at the specified
            % value of s.
            outVal = obj.conjugate(inS);
            outVal( inS > obj.limit ) = Inf;
        end
        
        function outDeriv = ConjugateDerivative( obj, inS )
            % ConjugateDerivative returns the derivative of the conjugate
            % at the specified value of s.
            outDeriv = obj.conjugateDerivative(inS);
            outDeriv( inS > obj.limit ) = NaN;
        end
        
        function outDeriv = SecondDerivativeAt1( obj )
            % SecondDerivative returns the second derivative of the
            % phi-divergence at the specified value of t.
            outDeriv = obj.phi2Derivative(1);
        end
        
        function SetComputationLimit( obj, inDistr, inRho )
            assert( all(inDistr >= 0) )
            if ~isinf( obj.limit )
                return
            end
            distr = inDistr/sum(inDistr);
            t = fsolve(@(t)obj.func(t) - inRho/min(distr(distr>0)), 2);
            s = fsolve(@(s)obj.conjugateDerivative(s) - t, 1);
            obj.computationLimit = s;
        end
        
        function rho = Rho( obj, alpha, obs )
            rho = obj.SecondDerivativeAt1() / (2*sum(obs)) * ...
                chi2inv(1-alpha,length(obs) - 1);
        end
        
        function alpha = Alpha( obj, rho, obs )
            alpha = 1-chi2cdf(2*sum(obs)/obj.SecondDerivativeAt1() * rho, ...
                length(obs)-1);
        end
    end
end

