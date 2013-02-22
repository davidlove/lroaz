% Solution stores a solution to an LRLP problem.  Solution store:
% * Value of all variables
% * thetas of the master and true problem
% * Dual values for the second stage problems
% * Whether mu was initially feasible
% * Whether the solution was in the trust region interior

classdef Solution < handle
    
    properties (GetAccess=public, SetAccess=private)
        solution
        lambda
        mu
        theta
        secondStageValues
        secondStageDuals
        muFeasible
        trustRegionInterior
    end
    
    properties (GetAccess=private, SetAccess=immutable)
        numVariables
        numScen
        numTheta
    end
    
    % Enumeration properties
    properties (GetAccess=private, SetAccess=immutable)
        MASTER
        TRUE
        SLOPE
        INTERCEPT
    end
    
    methods (Access=public)
        
        function obj = Solution( lp, cutType )
            if nargin < 2
                cutType = 'single';
            end
            if ~isa( lp, 'LPModel' )
                error('Solution must be initialized with an LPModel');
            end
            
            obj.MASTER = 1;
            obj.TRUE = 2;
            
            obj.SLOPE = 1;
            obj.INTERCEPT = 2;
            
            obj.numVariables = size(lp.A,2);            
            obj.numScen = lp.numScenarios;
            
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
            if numel(inX) ~= self.numVariables
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
        
        function SetTheta( self, inTheta, type )
            if numel(inTheta) ~= self.numScen
                error( 'Solution:SetTheta:size', ...
                    ['Theta has size ' num2str(numel(inTheta)) ...
                    ', should be ' num2str(self.numScen)] )
            end
            switch type
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
        
        function Reset( self )
            self.solution = zeros( self.numVariables, 1 );
            self.lambda = [];
            self.mu = [];
            
            self.theta = zeros( self.numTheta, 2 );
            
            self.secondStageValues = -Inf( self.numScen, 1 );
            self.secondStageDuals = cell( self.numScen, 2 );
            
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
        
    end
end