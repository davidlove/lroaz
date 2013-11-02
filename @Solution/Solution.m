classdef Solution < matlab.mixin.Copyable
    
    properties (GetAccess=private, SetAccess=private)
        solution
        lambda
        mu
        theta
        secondStageSolutions
        secondStageValues
        secondStageDuals
        muFeasible
        trustRegionInterior
    end
    
    properties (GetAccess=private, SetAccess=immutable)
        numVariables
        numScen
        numTheta
        phiLimit
        isObserved
    end
    
    % Enumeration properties
    properties (GetAccess=private, SetAccess=immutable)
        MASTER
        TRUE
        SLOPE
        INTERCEPT
    end
    
    methods (Access=public)
        
        function obj = Solution( varargin )
            % Solution stores a solution to an LRLP problem.  Solution store:
            % * Value of all variables
            % * thetas of the master and true problem
            % * Dual values for the second stage problems
            % * Whether mu was initially feasible
            % * Whether the solution was in the trust region interior
            %
            % Keyword arguments (*required arguments):
            %    lp*: LPModel object defining the SLP-2
            %    phi*: PhiDivergence object or string defining the phi-divergence
            %    cuttype: 'single' cut or 'multi' cut
            
            % -------------------------------------------------------------------
            % Input Parsing
            % -------------------------------------------------------------------
            
            requiredArgs = {'lp', 'obs', 'phi'};
            
            if mod(length(varargin),2) == 1
                error('Arguments must be key, value pairs')
            end
            
            for vv = 1:2:length(varargin)
                key = varargin{vv};
                value = varargin{vv+1};
                switch key
                    case 'cuttype'
                        cutType = value;
                    case 'lp'
                        lp = value;
                    case 'obs'
                        obs = value;
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
            
            if ~exist('cutType', 'var')
                cutType = 'multi';
            end
            
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
            
            if ~isa(lp, 'LPModel')
                error('lp must be an LPModel object')
            end
            
            if ~ismember( cutType, {'single','multi'} )
                error([cutType ' is not a valid cut type'])
            end
            
            % -------------------------------------------------------------------
            % -------------------------------------------------------------------
            
            obj.MASTER = 1;
            obj.TRUE = 2;
            
            obj.SLOPE = 1;
            obj.INTERCEPT = 2;
            
            obj.numVariables = size(lp.A,2);            
            obj.numScen = lp.numScenarios;
            obj.phiLimit = min( phi.limit, phi.computationLimit );
            obj.isObserved = obs(:) > 0;
            
            switch cutType
                case 'single'
                    obj.numTheta = 1;
                case 'multi'
                    obj.numTheta = obj.numScen;
                otherwise
                    error('Cut types available: single and multi')
            end
            
            obj.Reset();
        end
        
        function SetX( self, inX )
            if ~isempty(self.solution)
                error('Solution:SetX:setagain', 'X has already been set')
            elseif numel(inX) ~= self.numVariables
                error('Solution:SetX:size', ...
                    ['X has size ' num2str(numel(inX)) ...
                    ', should be ' num2str(self.numVariables)])
            end
            self.solution = inX(:);
        end
        
        function SetLambda( self, inLambda )
            if numel(inLambda) ~= 1
                error( 'Solution:SetLambda:size', ...
                    ['Lambda has size ' num2str(numel(inLambda)) ...
                    ', must have size 1'] )
            elseif inLambda < 0
                error( 'Solution:SetLambda:sign', ...
                    'Lambda must be non-negative' )
            end
            self.lambda = inLambda(:);
        end
        
        function SetMu( self, inMu )
           if numel(inMu) ~= 1
                error( 'Solution:SetMu:size', ...
                    ['Mu has size ' num2str(numel(inMu)) ...
                    ', must have size 1'] )
           end
            self.mu = inMu(:);
        end
        
        function SetTheta( self, inTheta, inType )
            if numel(inTheta) ~= size( self.theta, 1 );
                error( 'Solution:SetTheta:size', ...
                    ['Theta has size ' num2str(numel(inTheta)) ...
                    ', should be ' num2str(self.numScen)] )
            end
            switch inType
                case 'master'
                    typeN = self.MASTER;
                case 'true'
                    typeN = self.TRUE;
                otherwise
                    error('Solution:SetTheta:type', ...
                        'type must be ''master'' or ''true''')
            end
            self.theta(:,typeN) = inTheta(:);
        end
        
        function SetSecondStageValue( self, inScen, inValue )
            if inScen < 1 || inScen > self.numScen
                error( 'Solution:SetSecondStageValue:badscen', ...
                    ['Scenario number must be between 1 and ' ...
                    num2str(self.numScen)] )
            end
            self.secondStageValues(inScen) = inValue;
            if all(self.secondStageValues > -Inf)
                % NOTE: Strict inequality is required for most phi
                % divergences, but variationdistance allows for
                % <= phiLimit.
                tolerBound = 1e-6 * max(self.phiLimit, 1);
                self.muFeasible = all(self.S() <= self.phiLimit ...
                    + ~self.isObserved*tolerBound);
            end
        end
        
        function SetSecondStageDual( self, inScen, inDual, inType )
            if inScen < 1 || inScen > self.numScen
                error( 'Solution:SetSecondStageDual:badscen', ...
                    ['Scenario number must be between 1 and ' ...
                    num2str(self.numScen)] )
            end
            switch inType
                case 'slope'
                    type = self.SLOPE;
                    if numel(inDual) ~= self.numVariables
                        error( 'Solution:SetSecondStageDual:size', ...
                            ['Dual slope must have same number of ' ...
                            'variables as solution.'] )
                    end
                case 'int'
                    type = self.INTERCEPT;
                    if numel(inDual) ~= 1
                        error( 'Solution:SetSecondStageDual:size', ...
                        'Dual intercept must be a single element' )
                    end
                otherwise
                    error( 'Solution:SetSecondStageDual:type', ...
                        'Types: ''slope'' and ''int''' );
            end
            
            self.secondStageDuals{ inScen, type } = inDual(:).';
        end
        
        function SetSecondStageSolution( self, inScen, inSol)
            if inScen < 1 || inScen > self.numScen
                error( 'Solution:SetSecondStageValue:badscen', ...
                    ['Scenario number must be between 1 and ' ...
                    num2str(self.numScen)] )
            end
            self.secondStageSolutions{inScen} = inSol;
        end
            
        function SetTrustRegionInterior( self, inTF )
            if ~isa( inTF, 'logical' )
                error( 'Solution:SetTrustRegionInterior:logical', ...
                    'Must be a logical variable' )
            end
            if ~isempty( self.trustRegionInterior )
                error( 'Solution:SetTrustRegionInterior:setagain', ...
                    'Trust region status already set' )
            end
            self.trustRegionInterior = inTF;
        end
        
        function Reset( self )
            self.solution = [];
            self.lambda = [];
            self.mu = [];
            
            self.theta = -Inf( self.numTheta, 2 );
            
            self.secondStageValues = -Inf( self.numScen, 1 );
            self.secondStageDuals = cell( self.numScen, 2 );
            self.secondStageSolutions = cell( self.numScen, 1 );
            
            self.muFeasible = [];
            self.trustRegionInterior = [];
        end
        
        function outX = X( self )
            outX = self.solution;
        end
        
        function outL = Lambda( self )
            outL = self.lambda;
        end
        
        function outM = Mu( self )
            outM = self.mu;
        end
        
        function outT = ThetaMaster( self )
            outT = self.theta( :,self.MASTER );
        end
        
        function outT = ThetaTrue( self )
            outT = self.theta( :,self.TRUE );
        end
        
        function outS = S( self )
            if self.lambda ~= 0
                outS = (self.secondStageValues - self.mu) / self.lambda;
            else
                relDiff = (self.secondStageValues - self.mu) / abs(self.mu);
                outS = zeros(size(self.secondStageValues));
                tol = 1e-6;
                outS(relDiff < -tol) = -Inf;
                outS(relDiff > tol) = Inf;
            end
        end
        
        function outL = Limit( self )
            outL = self.phiLimit;
        end
        
        function outValues = SecondStageValues( self )
            outValues = self.secondStageValues;
        end
        
        function outSlope = SecondStageSlope( self, inScen )
            outSlope = self.secondStageDuals{ inScen, self.SLOPE };
        end
        
        function outInt = SecondStageIntercept( self, inScen )
            outInt = self.secondStageDuals{ inScen, self.INTERCEPT };
        end
        
        function outSol = SecondStageSolution( self, inScen )
            outSol = self.secondStageSolutions{ inScen };
        end
        
        function outTF = MuFeasible( self )
            outTF = self.muFeasible;
        end
        
        function outTF = TrustRegionInterior( self )
            outTF = self.trustRegionInterior;
        end
        
    end
end