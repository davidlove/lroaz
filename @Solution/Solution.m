% Solution stores a solution to an LRLP problem.  Solution store:
% * Value of all variables
% * thetas of the master and true problem
% * Dual values for the second stage problems
% * Whether mu was initially feasible
% * Whether the solution was in the trust region interior

classdef Solution < handle
    
    properties (GetAccess=public, SetAccess=public)
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
        
        function Reset( self )
            self.solution = zeros( self.numVariables, 1 );
            self.lambda = 0;
            self.mu = 0;
            
            self.theta = zeros( self.numTheta, 2 );
            
            self.secondStageValues = -Inf( self.numScen, 1 );
            self.secondStageDuals = cell( self.numScen, 2 );
            
            self.muFeasible = [];
            self.trustRegionInterior = [];
        end
        
        function xOut = X( self )
            xOut = self.solution;
        end
        
        function lOut = Lambda( self )
            lOut = self.lambda;
        end
        
        function mOut = Mu( self )
            mOut = self.mu;
        end
        
        function tOut = ThetaMaster( self )
            tOut = self.theta( :,self.MASTER );
        end
        
        function tOut = ThetaTrue( self )
            tOut = self.theta( :,self.TRUE );
        end
        
    end
end