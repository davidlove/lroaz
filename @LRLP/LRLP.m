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
        SLOPE
        INTERCEPT
    end
    
%     Bender's Decomposition Parameters
    properties (Access=private)
        candidateSolution
        bestSolution
        secondBestSolution
        thetaTrue
        zLower
        zUpper
        objectiveCutsMatrix
        objectiveCutsRHS
        feasibilityCutsMatrix
        feasibilityCutsRHS
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
        function obj = LRLP( inLPModel, inGammaPrime, inNumObsPerScen, inOptimizer )
            if nargin < 1 || ~isa(inLPModel,'LPModel')
                error('LRLP must be initialized with an LPModel as its first argument')
            end
            
            obj.lpModel = inLPModel;
            
            if nargin < 4
                inOptimizer = 'linprog';
                if nargin < 3
                    inNumObsPerScen = ones(1,obj.lpModel.numScenarios);
                    if nargin < 2
                        inGammaPrime = 0.5;
                    end
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
            obj.optimizer = inOptimizer;
            
            obj.LAMBDA = size( obj.lpModel.A, 2 ) + 1;
            obj.MU = obj.LAMBDA + 1;
            obj.THETA = obj.MU + 1;
            
            obj.SLOPE = 1;
            obj.INTERCEPT = 2;
            
            obj.InitializeBenders();
        end
        
        % InitializeBenders initializes all Bender's Decomposition
        % parameters 
        function InitializeBenders( obj )
            obj.objectiveCutsMatrix = [];
            obj.objectiveCutsRHS = [];
            obj.feasibilityCutsMatrix = [];
            obj.feasibilityCutsRHS = [];
            obj.zLower = -Inf;
            obj.zUpper = Inf;
            
            obj.ResetSecondStageSolutions();
            
            % Solve first stage LP
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
        function exitFlag = SolveMasterProblem( obj )
            obj.ResetSecondStageSolutions();
            
            cMaster = obj.GetMasterc();
            AMaster = obj.GetMasterA();
            bMaster = obj.GetMasterb();
            lMaster = obj.GetMasterl();
            uMaster = obj.GetMasteru();
            
            [obj.candidateSolution,~,exitFlag] = linprog( cMaster, ...
                [obj.objectiveCutsMatrix; obj.feasibilityCutsMatrix], ...
                [obj.objectiveCutsRHS   ; obj.feasibilityCutsRHS   ], ...
                AMaster, bMaster, ...
                lMaster, uMaster, ...
                [], [] );
            
            assert( length(obj.candidateSolution) == length(obj.bestSolution) );
        end
        
        % SolveSubProblems solves all subproblems, generates second stage
        % dual information
        function SolveSubProblems( obj )
            for scenarioNum = 1:obj.lpModel.numScenarios
                obj.SubProblem( scenarioNum );
            end
        end
        
        % GenerateCuts generates objective cut, and if necessary, generates
        % a feasibility cut and finds a good feasible value of mu
        function GenerateCuts( obj )
            if obj.Mu() <= max( obj.secondStageValues )
                obj.muIsFeasible = false;
                obj.GenerateFeasibilityCut();
                obj.FindFeasibleMu();
            else
                obj.muIsFeasible = true;
            end
            
            obj.GenerateObjectiveCut();
        end
        
        % UpdateTolerances updates the upper bound on the optimal value and
        % the objective and probability tolerances
        function UpdateTolerances( obj )
            
        end
        
    end
    
    methods (Access=private)
        % SubProblem solves an individual subproblem, updating the optimal
        % value and dual solution to the sub problem
        function SubProblem( obj, inScenNumber )
            q = obj.lpModel.Getq( inScenNumber );
            D = obj.lpModel.GetD( inScenNumber );
            d = obj.lpModel.Getd( inScenNumber );
            B = obj.lpModel.GetB( inScenNumber );
            l = obj.lpModel.Getl2( inScenNumber );
            u = obj.lpModel.Getu2( inScenNumber );
            [~,fval,~,~,pi] = linprog(q, ...
                [],[], ...
                D, d + B*obj.X(), ...
                l, u);
            
            obj.secondStageDuals{ inScenNumber, obj.SLOPE } = -pi.eqlin'*B;
            obj.secondStageDuals{ inScenNumber, obj.INTERCEPT } ...
                = - pi.eqlin'*d ...
                  - pi.upper(u<Inf)'*u(u<Inf) ...
                  - pi.lower(l~=0)'*l(l~=0);
            obj.secondStageValues( inScenNumber ) = fval;
        end
        
        % GenerateObjectiveCut generates an objective cut and adds it to
        % the matrix
        function GenerateObjectiveCut( obj )
            intermediateSlope = zeros(obj.lpModel.numScenarios, size(obj.lpModel.A,2)+2);
            
            for ii=1:obj.lpModel.numScenarios
                intermediateSlope(ii,:) ...
                    = [ obj.numObsTotal*obj.Lambda*(obj.secondStageDuals{ii,obj.SLOPE}./(obj.Mu - obj.secondStageValues(ii))), ...
                        obj.numObsTotal*log(obj.Lambda) + obj.numObsTotal - obj.numObsTotal*log(obj.Mu - obj.secondStageValues(ii)), ...
                       -obj.numObsTotal*obj.Lambda/(obj.Mu-obj.secondStageValues(ii))];
            end
            
            obj.thetaTrue = obj.GetExpectedSecondStage();
            
            slope = obj.numObsPerScen/obj.numObsTotal*intermediateSlope;
            intercept = obj.thetaTrue - slope*[obj.X;obj.Lambda;obj.Mu];
            
            obj.objectiveCutsMatrix = [obj.objectiveCutsMatrix; slope, -1];
            obj.objectiveCutsRHS = [obj.objectiveCutsRHS; -intercept];
        end
        
        % GenerateFeasibilityCut generates a feasibility cut and adds it to
        % the matrix
        function GenerateFeasibilityCut( obj )
            [~,hIndex] = max( obj.secondStageValues );
            
            feasSlope = [obj.secondStageDuals{hIndex,obj.SLOPE}, 0, -1, 0];
            feasInt = obj.secondStageDuals{hIndex,obj.INTERCEPT};
            
            obj.feasibilityCutsMatrix = [obj.feasibilityCutsMatrix; feasSlope];
            obj.feasibilityCutsRHS = [obj.feasibilityCutsRHS; -feasInt];
        end
        
        % FindFeasibleMu uses Newton's Method to find a feasible value of
        % mu
        function FindFeasibleMu( obj )
            [hMax hIndex] = max( obj.secondStageValues );
            mu = obj.lpModel.numScenarios/2 ...
                * obj.numObsPerScen(hIndex)*obj.Lambda() ...
                + hMax;
            
            for ii=1:200
                muOld = mu;
                mu = mu + (sum(obj.Lambda() * obj.numObsPerScen'./(mu - obj.secondStageValues))-1) / ...
                    sum(obj.Lambda() * obj.numObsPerScen'./((mu - obj.secondStageValues).^2));
                if abs(mu - muOld) < min(mu,muOld)*0.01
                    break
                end
            end
            obj.bestSolution( obj.MU ) = max(mu, hMax+1);
        end
        
        % GetExpectedSecondStage gets the expected value of the second
        % stage in the LRLP-2
        function outE = GetExpectedSecondStage( obj )
            outE = obj.numObsPerScen/obj.numObsTotal ...
                   * ( obj.numObsTotal*obj.Lambda*log(obj.Lambda) ...
                      - obj.numObsTotal*obj.Lambda*log(obj.Mu-obj.secondStageValues) );
        end
        
        % ResetSecondStageSolutions clears the second stage solution values
        % and dual solution information
        function ResetSecondStageSolutions( obj )
            obj.secondStageValues = -Inf(obj.lpModel.numScenarios,1);
            obj.secondStageDuals = cell(obj.lpModel.numScenarios,2);
            obj.candidateSolution = [];
            obj.thetaTrue = [];
        end
        
        % GetMasterc gets the cost vector for the master problem
        function cOut = GetMasterc( obj )
            cOut = [obj.lpModel.c, 0, 0, 0];
            cOut(obj.LAMBDA) = obj.nBar;
            cOut(obj.MU) = 1;
            cOut(obj.THETA) = 1;
        end
        
        % GetMasterA gets the constraint matrix for the master problem
        function AOut = GetMasterA( obj )
            AOut = [obj.lpModel.A, zeros( size(obj.lpModel.A,1), 3 )];
        end
        
        % GetMasterb gets the right hand side vector for the master problem
        function bOut = GetMasterb( obj )
            bOut = obj.lpModel.b;
        end
        
        % GetMasterl gets the lower bound for the master problem
        function lOut = GetMasterl( obj )
            lOut = [obj.lpModel.l; 0; 0; 0];
            lOut(obj.LAMBDA) = 0;
            lOut(obj.MU) = -Inf;
            lOut(obj.THETA) = -Inf;
        end
        
        % GetMasteru gets the upper bound for the master problem
        function uOut = GetMasteru( obj )
            uOut = [obj.lpModel.u; 0; 0; 0];
            uOut(obj.LAMBDA) = Inf;
            uOut(obj.MU) = Inf;
            uOut(obj.THETA) = Inf;
        end
        
    end
    
%     Accessor methods
    methods (Access=public)
        % X returns the best value of decisions x
        function outX = X( obj )
            outX = obj.bestSolution( 1:size(obj.lpModel.A,2) );
        end
        
        % Lambda returns the best value of lambda
        function outLambda = Lambda( obj )
            outLambda = obj.bestSolution( obj.LAMBDA );
        end
        
        % Mu returns the best value of mu
        function outMu = Mu( obj )
            outMu = obj.bestSolution( obj.MU );
        end
        
        % NumObjectiveCuts returns the number of objective cuts
        function outNum = NumObjectiveCuts( obj )
            outNum = size(obj.objectiveCutsMatrix,1);
        end
        
        % NumFeasibilityCuts returns the number of objective cuts
        function outNum = NumFeasibilityCuts( obj )
            outNum = size(obj.feasibilityCutsMatrix,1);
        end
        
    end
        
end