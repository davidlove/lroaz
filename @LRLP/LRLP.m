% LRLP solves a 2 Stage Likelihood Robust Linear Program with Recourse
% (LRLP-2) via the modified Bender's Decomposition proposed by David Love.
% This class uses the LPModel class to create and store the LP data.

classdef LRLP < handle
    
%     Problem Parameters
    properties (GetAccess=public, SetAccess=immutable)
        lpModel
        gammaPrime
        numObsPerScen
        numObsTotal
        nBar
        optimizer
    end
    
%     Enumeration Parameters
    properties (GetAccess=private, SetAccess=immutable)
        LAMBDA
        MU
        THETA
    end
    
%     Bender's Decomposition Parameters
    properties (Access=private)
        candidateSolution
        bestSolution
        secondBestSolution
        zLower
        zUpper
        objectiveCutsMatrix
        feasibilityCutsMatrix
        secondStageValues
        secondStageDuals
    end
    
%     LRLP Solution Parameters
    properties (GetAccess=public, SetAccess=private)
        pWorst
    end
    
%     Boolean parameters for the algorithm
    properties (Access=private)
        muIsFeasible
    end
    
    methods
        % LRLP Constructor checks initialization conditions and generates
        % all immutable properties for the algorithm
        function obj = LRLP( inLPModel, inGammaPrime, inNumObsPerScen )
            if nargin < 1 || ~isa(inLPModel,'LPModel')
                error('LRLP must be initialized with an LPModel as its first argument')
            end
            
            obj.lpModel = inLPModel;
            
            if nargin < 3
                inNumObsPerScen = ones(1,obj.lpModel.numScenarios);
                if nargin < 2
                    inGammaPrime = 0.5;
                end
            end
            
            if obj.lpModel.numStages ~= 2
                error('Must use a 2 Stage LP')
            end
            if length(inNumObsPerScen) ~= obj.lpModel.numScenarios
                error('Size of observations differs from number of scenarios')
            end
            if inGammaPrime < 0 || inGammaPrime > 1
                error('Gamma Prime must be between 0 and 1')
            end
            
            obj.gammaPrime = inGammaPrime;
            obj.numObsPerScen = inNumObsPerScen;
            obj.numObsTotal = sum(obj.numObsPerScen);
            obj.nBar = obj.numObsTotal*(log(obj.numObsTotal)-1) ...
                - log(obj.gammaPrime);
            obj.optimizer = 'linprog';
            
            obj.LAMBDA = size( obj.lpModel.A, 2 ) + 1;
            obj.MU = obj.LAMBDA + 1;
            obj.THETA = obj.MU + 1;
            
            obj.InitializeBenders();
        end
        
        % InitializeBenders initializes all Bender's Decomposition
        % parameters 
        function InitializeBenders( obj )
            obj.objectiveCutsMatrix = [];
            obj.feasibilityCutsMatrix = [];
            obj.secondStageValues = cell(1,obj.lpModel.numScenarios);
            obj.secondStageDuals = cell(1,obj.lpModel.numScenarios);
            obj.zLower = -Inf;
            obj.zUpper = Inf;
            
            % Solve first stage LP
            obj.candidateSolution = [];
            switch obj.optimizer
                case 'linprog'
                    [x0,~,exitFlag] = linprog(obj.lpModel.c, ...
                        [], [], ...
                        obj.lpModel.A, obj.lpModel.b, ...
                        obj.lpModel.l, obj.lpModel.u);
                otherwise
                    error(['Optimizer ' obj.optimizer ' is not defined'])
            end
            if exitFlag ~= 1
                error('Could not solve first stage LP')
            end
            
            obj.bestSolution = [x0; 0; 0; 0];
            obj.bestSolution(obj.LAMBDA) = 1;
            obj.bestSolution(obj.MU) = -Inf;
            obj.bestSolution(obj.THETA) = -Inf;
            
            assert( length(obj.bestSolution) == obj.THETA );
            
            obj.SolveSubProblems()
            obj.GenerateCuts()
            obj.UpdateTolerances()
        end
            
                
        % SolveMasterProblem clears candidate solution and second stage
        % information, then solves the master problem
        function SolveMasterProblem( obj )
            
        end
        
        % SolveSubProblems solves all subproblems, generates second stage
        % dual information
        function SolveSubProblems( obj )
            
        end
        
        % GenerateCuts generates objective cut, and if necessary, generates
        % a feasibility cut and updates the value of mu
        function GenerateCuts( obj )
            
        end
        
        % UpdateTolerances updates the upper bound on the optimal value and
        % the objective and probability tolerances
        function UpdateTolerances( obj )
            
        end
        
    end
    
    methods (Access=private)
        % H solves an individual subproblem
        function H( obj, inScenNumber )
            
        end
        
        % GenerateObjectiveCut generates an objective cut and adds it to
        % the matrix
        function GenerateObjectiveCut( obj )
            
        end
        
        % GenerateFeasibilityCut generates a feasibility cut and adds it to
        % the matrix
        function GenerateFeasibilityCut( obj )
            
        end
        
        % FindFeasibleMu uses Newton's Method to find a feasible value of
        % mu
        function FindFeasibleMu( obj )
            
        end
        
    end
        
end