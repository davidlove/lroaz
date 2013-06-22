classdef PhiDivergence
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
                    obj.func = @(t) t.*log(t) - t + 1;
                    obj.conjugate = @(s) exp(s) - 1;
                    obj.conjugateDerivative = @(s) exp(s);
                    obj.phi2Derivative = @(t) 1./t;
                    obj.limit = Inf;
                case 'chi2'
                    obj.func = @(t) (t-1).^2./t;
                    obj.conjugate = @(s) 2 - 2*sqrt(1-s);
                    obj.conjugateDerivative = @(s) 1./sqrt(1-s);
                    obj.phi2Derivative = @(t) 2./t^3;
                    obj.limit = 1;
                case 'mchi2'
                    obj.func = @(t) (t-1).^2;
                    obj.conjugate = @(s) max(-1, s + s.^2/4);
                    obj.conjugateDerivative = @(s) (1+s/2).*(s >= -2);
                    obj.phi2Derivative = @(t) 2;
                    obj.limit = Inf;
                case 'hellinger'
                    obj.func = @(t) (sqrt(t)-1).^2;
                    obj.conjugate = @(s) s./(1-s);
                    obj.conjugateDerivative = @(s) 1./((s-1).^2);
                    obj.phi2Derivative = @(t) 1./(2*t^(3/2));
                    obj.limit = 1;
                otherwise
                    error('PhiDivergence:PhiDivergence:Unknown', ...
                        ['Unknown phi type ' inDivergence])
            end
            obj.divergence = inDivergence;
        end
        
        function outVal = Value( obj, inT )
            outVal = obj.func(inT);
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
    end
end

