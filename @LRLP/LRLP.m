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
    properties %(Access=private)
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
    properties %(Access=private)
        candidateMuIsFeasible
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
            
            obj.candidateSolution = [x0; 0; 0; 0];
            obj.candidateSolution(obj.LAMBDA) = 1;
            obj.candidateSolution(obj.MU) = -Inf;
            obj.candidateSolution(obj.THETA) = -Inf;
            
            assert( length(obj.candidateSolution) == obj.THETA );
            
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
        % dual information, determines whether mu is feasible
        function SolveSubProblems( obj, solution )
            if nargin < 2
                solution = obj.candidateSolution;
            end
            
            for scenarioNum = 1:obj.lpModel.numScenarios
                obj.SubProblem( scenarioNum, solution );
            end
            
            if obj.GetMu( solution ) > max( obj.secondStageValues )
                obj.candidateMuIsFeasible = true;
            else
                obj.candidateMuIsFeasible = false;
            end
        end
        
        % GenerateCuts generates objective cut, and if necessary, generates
        % a feasibility cut and finds a good feasible value of mu
        function GenerateCuts( obj )
            if ~obj.candidateMuIsFeasible
                obj.GenerateFeasibilityCut();
                obj.FindFeasibleMu();
            end
            
            obj.FindExpectedSecondStage();
            
            obj.GenerateObjectiveCut();
        end
        
        % UpdateTolerances updates the upper bound on the optimal value and
        % the objective and probability tolerances
        function UpdateTolerances( obj )
            
        end
        
        % Plot produces plots of a decision variable, lambda and mu all
        % together
        function Plot( obj, variableNumber )
            if nargin < 2
                variableNumber = 1;
            end
            
            if variableNumber > size(obj.lpModel.A,2)
                error('Must choose variable in the original LP model')
            end
            
            origCandidate = obj.candidateSolution;
            origBest = obj.bestSolution;
            origSecondBest = obj.secondBestSolution;
            origSecondDuals = obj.secondStageDuals;
            origSecondValues = obj.secondStageValues;
            origTheta = obj.thetaTrue;
            
            obj.PlotStep( variableNumber );
            obj.PlotStep( obj.LAMBDA );
            obj.PlotStep( obj.MU );
            
            % Ensure that no second stage properties while plotting
            assert( isequal( origCandidate, obj.candidateSolution ) )
            assert( isequal( origBest, obj.bestSolution ) )
            assert( isequal( origSecondBest, obj.secondBestSolution ) )
            assert( isequal( origSecondDuals, obj.secondStageDuals ) )
            assert( isequal( origSecondValues, obj.secondStageValues ) )
            assert( isequal( origTheta, obj.thetaTrue ) )
            
        end
        
    end
    
    methods (Access=private)
        
        PlotStep( obj, inVariableNumber );
        PlotFeasibilityCut( obj, inVariableNumber, inCutNumber, ylim );
        PlotObjectiveCut( obj, inVariableNumber, inCutNumber, inBounds )
        
        % SubProblem solves an individual subproblem, updating the optimal
        % value and dual solution to the sub problem
        function SubProblem( obj, inScenNumber, inSolution )
            q = obj.lpModel.Getq( inScenNumber );
            D = obj.lpModel.GetD( inScenNumber );
            d = obj.lpModel.Getd( inScenNumber );
            B = obj.lpModel.GetB( inScenNumber );
            l = obj.lpModel.Getl2( inScenNumber );
            u = obj.lpModel.Getu2( inScenNumber );
            
            xLocal = obj.GetX( inSolution );
            
            [~,fval,~,~,pi] = linprog(q, ...
                [],[], ...
                D, d + B*xLocal, ...
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
            xLocal = obj.GetX( obj.candidateSolution );
            lambdaLocal = obj.GetLambda( obj.candidateSolution );
            muLocal = obj.GetMu( obj.candidateSolution );
            
            intermediateSlope = zeros(obj.lpModel.numScenarios, size(obj.lpModel.A,2)+2);
            
            for ii=1:obj.lpModel.numScenarios
                intermediateSlope(ii,:) ...
                    = [ obj.numObsTotal*lambdaLocal*(obj.secondStageDuals{ii,obj.SLOPE}./(muLocal - obj.secondStageValues(ii))), ...
                        obj.numObsTotal*log(lambdaLocal) + obj.numObsTotal - obj.numObsTotal*log(muLocal - obj.secondStageValues(ii)), ...
                       -obj.numObsTotal*lambdaLocal/(muLocal-obj.secondStageValues(ii))];
            end
            
            slope = obj.numObsPerScen/obj.numObsTotal*intermediateSlope;
            intercept = obj.thetaTrue - slope*[xLocal;lambdaLocal;muLocal];
            
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
            lambdaLocal = obj.GetLambda( obj.candidateSolution );
            
            [hMax hIndex] = max( obj.secondStageValues );
            mu = obj.lpModel.numScenarios/2 ...
                * obj.numObsPerScen(hIndex)*lambdaLocal() ...
                + hMax;
            
            for ii=1:200
                muOld = mu;
                mu = mu + (sum(lambdaLocal() * obj.numObsPerScen'./(mu - obj.secondStageValues))-1) / ...
                    sum(lambdaLocal() * obj.numObsPerScen'./((mu - obj.secondStageValues).^2));
                if abs(mu - muOld) < min(mu,muOld)*0.01
                    break
                end
            end
            obj.candidateSolution( obj.MU ) = max(mu, hMax+1);
        end
        
        % FindExpectedSecondStage gets the expected value of the second
        % stage in the LRLP-2
        function FindExpectedSecondStage( obj, inSolution )
            if nargin < 2
                inSolution = obj.candidateSolution;
            end
            
            lambdaLocal = obj.GetLambda( inSolution );
            muLocal = obj.GetMu( inSolution );
            
            obj.thetaTrue = obj.numObsPerScen/obj.numObsTotal ...
                   * ( obj.numObsTotal*lambdaLocal*log(lambdaLocal) ...
                      - obj.numObsTotal*lambdaLocal*log(muLocal-obj.secondStageValues) );
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
        
        % GetX gets the decisions x from the given solution
        function outX = GetX( obj, solution )
            outX = solution( 1:size(obj.lpModel.A,2) );
        end
        
        % GetLambda gets lambda from the given solution
        function outLambda = GetLambda( obj, solution )
            outLambda = solution( obj.LAMBDA );
        end
        
        % GetMU gets mu from the given solution
        function outMu = GetMu( obj, solution )
            outMu = solution( obj.MU );
        end
        
    end
    
%     Accessor methods
    methods (Access=public)
        % X returns the best value of decisions x
        function outX = X( obj )
            outX = obj.GetX( obj.bestSolution );
        end
        
        % Lambda returns the best value of lambda
        function outLambda = Lambda( obj )
            outLambda = obj.GetLambda( obj.bestSolution );
        end
        
        % Mu returns the best value of mu
        function outMu = Mu( obj )
            outMu = obj.GetMu( obj.bestSolution );
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
