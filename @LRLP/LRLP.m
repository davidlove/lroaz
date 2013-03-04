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
        objectiveTolerance
        probabilityTolerance
        trustRegionRhoBound
        trustRegionScaleDown
        trustRegionScaleUp
        trustRegionEta
        trustRegionRatio
        trustRegionMinSize
        trustRegionMaxSize
    end
    
    %     Enumeration Parameters
    properties (GetAccess=private, SetAccess=immutable)
        LAMBDA
        MU
        THETA
        SLOPE
        INTERCEPT
        NO_SCALE
        SCALE_UP
        SCALE_DOWN
    end
    
    %     Bender's Decomposition Parameters
    properties (Access=private)
        candidateSolution
        secondBestSolution
        zLower
        zUpper
        objectiveCutsMatrix
        objectiveCutsRHS
        feasibilityCutsMatrix
        feasibilityCutsRHS
        trustRegionSize
        trustRegionLower
        trustRegionUpper
        trustRegionRho
        objectiveScale
        optimizerOptions
    end
    
    %     LRLP Solution Parameters
    properties (GetAccess=public, SetAccess=private)
        bestSolution
        currentObjectiveTolerance
        currentProbabilityTolerance
        pWorst
        relativeLikelihood
    end
    
    %     Boolean parameters for the algorithm
    properties %(Access=private)
        newSolutionAccepted
        zLowerUpdated
        trustRegionScaled
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
            
            obj.objectiveTolerance = 1e-5;
            obj.probabilityTolerance = 1e-3;
            
            obj.trustRegionRhoBound = 1/4;
            obj.trustRegionScaleDown = 1/4;
            obj.trustRegionScaleUp = 3;
            obj.trustRegionEta = 1/50;
            obj.trustRegionRatio = 0.99;
            
            assert( obj.trustRegionEta <= obj.trustRegionRhoBound )
            
            obj.LAMBDA = size( obj.lpModel.A, 2 ) + 1;
            obj.MU = obj.LAMBDA + 1;
            obj.THETA = obj.MU + 1;
            
            obj.SLOPE = 1;
            obj.INTERCEPT = 2;
            
            obj.NO_SCALE = 0;
            obj.SCALE_UP = 1;
            obj.SCALE_DOWN = 2;
            
            obj.candidateSolution = Solution( obj.lpModel, 'single' );
            obj.InitializeBenders();
            
            obj.trustRegionMinSize = obj.trustRegionSize / 10;
            obj.trustRegionMaxSize = obj.trustRegionSize * Inf;%1000;
            assert( obj.trustRegionMinSize <= obj.trustRegionMaxSize );
        end
        
        % InitializeBenders initializes all Bender's Decomposition
        % parameters
        function InitializeBenders( obj )
            switch obj.optimizer
                case 'linprog'
                    obj.optimizerOptions = optimset( obj.optimizer );
                    obj.optimizerOptions = optimset( obj.optimizerOptions, ...
                        'MaxIter', 85 );
                case 'cplexlp'
                    % No options set yet for cplexlp
                otherwise
                    error(['Optimizer ' obj.optimizer ' is not defined'])
            end
            
            obj.objectiveCutsMatrix = [];
            obj.objectiveCutsRHS = [];
            obj.feasibilityCutsMatrix = [];
            obj.feasibilityCutsRHS = [];
            obj.zLower = -Inf;
            obj.zUpper = Inf;
            
            obj.bestSolution = [];
            obj.secondBestSolution = [];
            
            obj.ResetSecondStageSolutions();
            obj.newSolutionAccepted = true;
            
            % Solve first stage LP
            switch obj.optimizer
                case 'linprog'
                    [x0,~,exitFlag] = linprog(obj.lpModel.c, ...
                        [], [], ...
                        obj.lpModel.A, obj.lpModel.b, ...
                        obj.lpModel.l, obj.lpModel.u);
                case 'cplexlp'
                    [x0,~,exitFlag] = cplexlp(obj.lpModel.c, ...
                        [], [], ...
                        obj.lpModel.A, obj.lpModel.b, ...
                        obj.lpModel.l, obj.lpModel.u);
                otherwise
                    error(['Optimizer ' obj.optimizer ' is not defined'])
            end
            if exitFlag ~= 1
                error('Could not solve first stage LP')
            end
            
            obj.objectiveScale = 1;
%             obj.objectiveScale = obj.nBar ...
%                 / median(obj.lpModel.c(obj.lpModel.c ~= 0));
            
            obj.trustRegionSize = 1*max(abs(x0));
            
            obj.candidateSolution.SetX( x0 );
            obj.candidateSolution.SetLambda( 1 );
            obj.candidateSolution.SetMu( -Inf );
            
            obj.SolveSubProblems();
            obj.GenerateCuts();
            
            obj.UpdateBestSolution();
            obj.UpdateTolerances();
        end
        
        
        % SolveMasterProblem clears candidate solution and second stage
        % information, solves the master problem and determines whether the
        % solution is in the interior of the trust region
        function exitFlag = SolveMasterProblem( obj )
            obj.ResetSecondStageSolutions();
            
            cMaster = obj.GetMasterc();
            AMaster = obj.GetMasterA();
            bMaster = obj.GetMasterb();
            lMaster = obj.GetMasterl();
            uMaster = obj.GetMasteru();
            
            currentBest = obj.GetDecisions( obj.bestSolution );

            % Best solution should be feasible
            feasTol = max(abs(currentBest)) * 1e-8;
            assert( all( abs( AMaster*currentBest - bMaster ) ...
                < feasTol ), ...
                ['bestSolution master infeasibility: ' ...
                num2str(max(abs( AMaster*currentBest - bMaster )))])
            assert( all( obj.objectiveCutsMatrix*currentBest ...
                - obj.objectiveCutsRHS < feasTol ), ...
                ['bestSolution objective infeasibility: ' ...
                num2str( max(obj.objectiveCutsMatrix*currentBest ...
                - obj.objectiveCutsRHS) )] )
            assert( all( obj.feasibilityCutsMatrix*currentBest ...
                - obj.feasibilityCutsRHS < feasTol ), ...
                ['bestSolution feasibility infeasibility: ' ...
                num2str( max(obj.feasibilityCutsMatrix*currentBest ...
                - obj.feasibilityCutsRHS) )] )

            obj.candidateSolution.Reset;
            
            switch obj.optimizer
                case 'linprog'
                    [currentCandidate,~,exitFlag] = linprog( cMaster, ...
                        [obj.objectiveCutsMatrix; obj.feasibilityCutsMatrix], ...
                        [obj.objectiveCutsRHS   ; obj.feasibilityCutsRHS   ], ...
                        AMaster, bMaster, ...
                        lMaster, uMaster, ...
                        [], obj.optimizerOptions );
                case 'cplexlp'
                    [currentCandidate,~,exitFlag] = cplexlp( cMaster, ...
                        [obj.objectiveCutsMatrix; obj.feasibilityCutsMatrix], ...
                        [obj.objectiveCutsRHS   ; obj.feasibilityCutsRHS   ], ...
                        AMaster, bMaster, ...
                        lMaster, uMaster, ...
                        [], obj.optimizerOptions );
                otherwise
                    error(['Optimizer ' obj.optimizer ' is not defined'])
            end

            obj.candidateSolution.SetX( currentCandidate(1:end-3) )
            obj.candidateSolution.SetLambda( currentCandidate(obj.LAMBDA) )
            obj.candidateSolution.SetMu( currentCandidate(obj.MU) )
            obj.candidateSolution.SetTheta( currentCandidate(obj.THETA), 'master' )
            
            % If Matlab fails to find an optimal solution, whether it
            % recongnizes it (bet setting exitFlag < 1) or  not (by not
            % beating bestSolution), return immediately and let the program
            % controlling LRLP handle it.
            if exitFlag ~= 1 || ...
                    cMaster*(currentBest - currentCandidate) < 0
                if exitFlag == 1
                    exitFlag = -50;
                end
                return
            end
            
            % Any accepted solution should be better than the previous
            % best.
            assert( exitFlag ~= 1 || ...
                cMaster * (currentBest - currentCandidate) >= 0, ...
                ['Actual objective drop = ' ...
                num2str( cMaster * (currentBest - currentCandidate) )])
            
            upperT = ((1+obj.trustRegionRatio)*obj.trustRegionUpper ...
                + (1-obj.trustRegionRatio)*obj.trustRegionLower) ...
                / 2;
            lowerT = ((1-obj.trustRegionRatio)*obj.trustRegionUpper ...
                + (1+obj.trustRegionRatio)*obj.trustRegionLower) ...
                / 2;
            
            ind = 1:obj.THETA-1;
            
            assert( all( upperT(ind) - lowerT(ind) ...
                - obj.trustRegionRatio*( obj.trustRegionUpper(ind) - obj.trustRegionLower(ind) ) ...
                < 1e-6 ) )
            assert( all( (upperT(ind) + lowerT(ind))/2 ...
                - (obj.trustRegionUpper(ind) + obj.trustRegionLower(ind))/2 ...
                < 1e-6 ) )
            
            obj.candidateSolution.SetTrustRegionInterior( ~any( lowerT(ind) > currentCandidate(ind) ...
                | currentCandidate(ind) > upperT(ind) ) );
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
        end
        
        % GenerateCuts generates objective cut, and if necessary, generates
        % a feasibility cut and finds a good feasible value of mu
        function GenerateCuts( obj )
            if ~obj.candidateSolution.MuFeasible
                obj.GenerateFeasibilityCut();
                obj.FindFeasibleMu();
            end
            
            obj.FindExpectedSecondStage();
            
            obj.GenerateObjectiveCut();
        end
        
        % UpdateTrustRegionSize updates the size of the trust region
        function UpdateTrustRegionSize( obj )
            obj.CalculateTrustRegionDecisions();
            
            switch obj.trustRegionScaled
                case obj.SCALE_DOWN
                    obj.trustRegionSize = max( obj.trustRegionMinSize, ...
                        obj.trustRegionScaleDown * obj.trustRegionSize );
                case obj.SCALE_UP
                    obj.trustRegionSize = min( obj.trustRegionMaxSize, ...
                        obj.trustRegionScaleUp * obj.trustRegionSize );
                case obj.NO_SCALE
                    % Do nothing
                otherwise
                    error(['Unknown trust region scaling result: ' ...
                        num2str( obj.trustRegionScaled )])
            end
        end
        
        % UpdateSolutions updates the best solutions, and the upper and
        % lower bounds
        function UpdateSolutions( obj )
            obj.CalculateTrustRegionDecisions();
            
            cMaster = obj.GetMasterc();
            if obj.newSolutionAccepted
                obj.UpdateBestSolution();
                bestDecisions = obj.GetDecisions( obj.bestSolution, 'true' );
                assert( cMaster*bestDecisions <= obj.zUpper, ...
                    ['rho = ' num2str(obj.trustRegionRho)])
                obj.zUpper = cMaster*bestDecisions;
            end
            
            if obj.newSolutionAccepted ...
                    && obj.candidateSolution.MuFeasible ...
                    && obj.candidateSolution.TrustRegionInterior
                candidateDecisions = obj.GetDecisions( obj.candidateSolution, 'master' );
                assert( cMaster*candidateDecisions >= obj.zLower )
                obj.zLowerUpdated = true;
                obj.zLower = cMaster*candidateDecisions;
            else
                obj.zLowerUpdated = false;
            end
        end
        
        % UpdateTolerances updates the upper bound on the optimal value and
        % the objective and probability tolerances
        function UpdateTolerances( obj )
            % obj.CalculateProbability();
            
            if obj.zLower > -Inf
                obj.currentObjectiveTolerance = (obj.zUpper - obj.zLower) ...
                    / min(abs(obj.zUpper),abs(obj.zLower));
            else
                obj.currentObjectiveTolerance = Inf;
            end
            obj.currentProbabilityTolerance = abs(1-sum(obj.pWorst));
        end
        
        function WriteProgress( obj )
            disp(' ')
            disp(['gamma'' = ' num2str(obj.gammaPrime)])
            disp(['Observations: ' num2str(obj.numObsPerScen)])
            disp([num2str(obj.NumObjectiveCuts) ' objective cuts, '...
                num2str(obj.NumFeasibilityCuts) ' feasibility cuts.'])
            
            if obj.candidateSolution.MuFeasible
                disp('No feasibility cut generated')
            else
                disp('Feasibility cut generated')
            end
            
            if obj.candidateSolution.TrustRegionInterior
                disp(['Candidate solution in trust region interior, '...
                    'rho = ' num2str(obj.trustRegionRho)])
            else
                disp(['Candidate solution on trust region boundary, '...
                    'rho = ' num2str(obj.trustRegionRho)])
            end
            
            switch obj.trustRegionScaled
                case obj.NO_SCALE
                    disp('No change in trust region scale')
                case obj.SCALE_UP
                    disp(['Trust region scaled up to ' ...
                        num2str(obj.trustRegionSize)])
                case obj.SCALE_DOWN
                    disp(['Trust region scaled down to ' ...
                        num2str(obj.trustRegionSize)])
                otherwise
                    error(['Unknown trust region result ' ...
                        num2str(obj.trustRegionScaled)])
            end
            
            if obj.newSolutionAccepted
                disp(['New solution, zupper = ' num2str(obj.zUpper)])
            else
                disp('No new solution accepted')
            end
            
            if obj.zLowerUpdated
                disp(['New lower bound, zlower = ' num2str(obj.zLower)])
            else
                disp('No new lower bound found')
            end
            
            disp(['Objective tolerance ' num2str(obj.currentObjectiveTolerance)])
            disp(['Probability tolerance ' num2str(obj.currentProbabilityTolerance)])
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
            
            obj.PlotStep( variableNumber );
            obj.PlotStep( obj.LAMBDA );
            obj.PlotStep( obj.MU );
        end
        
        % DoubleIterations doubles the maximum number of iterations
        % available to the linear programming solver
        function DoubleIterations( obj )
            switch obj.optimizer
                case 'linprog'
                    newIter = min( 1000, 2*optimget(obj.optimizerOptions,'MaxIter'));
                    obj.optimizerOptions = optimset( obj.optimizerOptions, ...
                        'MaxIter', newIter );
                case 'cplexlp'
                    
                otherwise
                    error(['Unknown optimizer ' obj.optimizer])
            end
            disp(['Max Iterations = ' num2str(newIter)])
        end
        
        % DeleteOldestCut deletes the oldest objective cut from the matrix
        % and right hand side vector
        function DeleteOldestCut( obj )
            obj.objectiveCutsMatrix = obj.objectiveCutsMatrix(2:end,:);
            obj.objectiveCutsRHS = obj.objectiveCutsRHS(2:end);
            assert( ~isempty(obj.objectiveCutsRHS) )
        end
        
        function DeleteOldestFeasibilityCut( obj )
            obj.feasibilityCutsMatrix = obj.feasibilityCutsMatrix(2:end,:);
            obj.feasibilityCutsRHS = obj.feasibilityCutsRHS(2:end);
            assert( ~isempty(obj.feasibilityCutsRHS) )
        end
        
    end
    
    methods (Access=private)
        
        PlotStep( obj, inVariableNumber );
        PlotFeasibilityCut( obj, inVariableNumber, inCutNumber, ylim );
        PlotObjectiveCut( obj, inVariableNumber, inCutNumber, inBounds )
        
        % SubProblem solves an individual subproblem, updating the optimal
        % value and dual solution to the sub problem
        function SubProblem( obj, inScenNumber, inSolution )
            q = obj.lpModel.Getq( inScenNumber ) * obj.objectiveScale;
            D = obj.lpModel.GetD( inScenNumber );
            d = obj.lpModel.Getd( inScenNumber );
            B = obj.lpModel.GetB( inScenNumber );
            l = obj.lpModel.Getl2( inScenNumber );
            u = obj.lpModel.Getu2( inScenNumber );
            
            xLocal = inSolution.X;
            
            options = optimset('Display','off');
            
            switch obj.optimizer
                case 'linprog'
                    [~,fval,~,~,pi] = linprog(q, ...
                        [],[], ...
                        D, d + B*xLocal, ...
                        l, u, ...
                        [], options);
                case 'cplexlp'
                    [~,fval,~,~,pi] = cplexlp(q, ...
                        [],[], ...
                        D, d + B*xLocal, ...
                        l, u, ...
                        [], []);
                otherwise
                    error(['Optimizer ' obj.optimizer ' is not defined'])
            end

            inSolution.SetSecondStageDual( inScenNumber, sparse(-pi.eqlin'*B), 'slope' )
            inSolution.SetSecondStageDual( inScenNumber, - pi.eqlin'*d ...
                - pi.upper(u<Inf)'*u(u<Inf) ...
                - pi.lower(l~=0)'*l(l~=0), 'int' )
            inSolution.SetSecondStageValue( inScenNumber, fval )
        end
        
        % GenerateObjectiveCut generates an objective cut and adds it to
        % the matrix
        function GenerateObjectiveCut( obj )
            xLocal = obj.candidateSolution.X;
            lambdaLocal = obj.candidateSolution.Lambda;
            muLocal = obj.candidateSolution.Mu;

            localValues = obj.candidateSolution.SecondStageValues;
            
            intermediateSlope = zeros(obj.lpModel.numScenarios, size(obj.lpModel.A,2)+2);
            
            for ii=1:obj.lpModel.numScenarios
                intermediateSlope(ii,:) ...
                    = [ obj.numObsTotal*lambdaLocal*(obj.candidateSolution.SecondStageSlope(ii)./(muLocal - localValues(ii))), ...
                    obj.numObsTotal*log(lambdaLocal) + obj.numObsTotal - obj.numObsTotal*log(muLocal - localValues(ii)), ...
                    -obj.numObsTotal*lambdaLocal/(muLocal-localValues(ii))];
            end
            
            slope = obj.numObsPerScen/obj.numObsTotal*intermediateSlope;
            intercept = obj.candidateSolution.ThetaTrue - slope*[xLocal;lambdaLocal;muLocal];
            
            obj.objectiveCutsMatrix = [obj.objectiveCutsMatrix; sparse([slope, -1])];
            obj.objectiveCutsRHS = [obj.objectiveCutsRHS; -intercept];
        end
        
        % GenerateFeasibilityCut generates a feasibility cut and adds it to
        % the matrix
        function GenerateFeasibilityCut( obj )
            [~,hIndex] = max( obj.candidateSolution.SecondStageValues );
            
            feasSlope = [obj.candidateSolution.SecondStageSlope(hIndex), 0, -1, 0];
            feasInt = obj.candidateSolution.SecondStageIntercept(hIndex);
            
            obj.feasibilityCutsMatrix = [obj.feasibilityCutsMatrix; sparse(feasSlope)];
            obj.feasibilityCutsRHS = [obj.feasibilityCutsRHS; -feasInt];
        end
        
        % FindFeasibleMu uses Newton's Method to find a feasible value of
        % mu
        function FindFeasibleMu( obj )
            lambdaLocal = obj.candidateSolution.Lambda;
            
            localValues = obj.candidateSolution.SecondStageValues;
            [hMax hIndex] = max( localValues );
            mu = obj.lpModel.numScenarios/2 ...
                * obj.numObsPerScen(hIndex)*lambdaLocal() ...
                + hMax;
            
            for ii=1:200
                muOld = mu;
                mu = mu + (sum(lambdaLocal() * obj.numObsPerScen'./(mu - localValues))-1) / ...
                    sum(lambdaLocal() * obj.numObsPerScen'./((mu - localValues).^2));
                if abs(mu - muOld) < min(mu,muOld)*0.01
                    break
                end
            end
            
            if mu > hMax
                obj.candidateSolution.SetMu( mu );
            else
                obj.candidateSolution.SetMu( hMax*1.001 );
            end
        end
        
        % FindExpectedSecondStage gets the expected value of the second
        % stage in the LRLP-2
        function FindExpectedSecondStage( obj, inSolution )
            if nargin < 2
                inSolution = obj.candidateSolution;
            end
            
            if isempty( inSolution.MuFeasible )
                error(['Must determine whether candidate mu is feasible ' ...
                    'before finding expected second stage value'])
            end
            if ~all( inSolution.SecondStageValues > -Inf )
                error('Must set second stage values before calculating expectation')
            end
            
            lambdaLocal = inSolution.Lambda;
            muLocal = inSolution.Mu;
            
            inSolution.SetTheta( obj.numObsPerScen/obj.numObsTotal ...
                * ( obj.numObsTotal*lambdaLocal*log(lambdaLocal) ...
                - obj.numObsTotal*lambdaLocal*log(muLocal-inSolution.SecondStageValues) ), ...
                'true' );
        end
        
        % CalculateRho calculates the value of trustRegionRho
        function CalculateTrustRegionDecisions( obj )
            cMaster = obj.GetMasterc();
            
            initialSolution = obj.GetDecisions( obj.bestSolution, 'true' );
            candidate = obj.GetDecisions( obj.candidateSolution, 'master' );
            truth = obj.GetDecisions( obj.candidateSolution, 'true' );
            
            trueDrop = cMaster*(initialSolution - truth);
            predictedDrop = cMaster*(initialSolution - candidate);
            
            % Infeasible candidate solutions might not produce a predicted
            % drop after they are made feasible
            assert( ~obj.candidateSolution.MuFeasible || predictedDrop >= 0 )
            
            if obj.candidateSolution.MuFeasible
                if isempty( obj.candidateSolution.TrustRegionInterior )
                    error( ['Undetermined whether solution exists in trust ' ...
                        'region interior'])
                end
                
                if trueDrop < obj.trustRegionRhoBound * predictedDrop
                    obj.trustRegionScaled = obj.SCALE_DOWN;
                elseif trueDrop > (1-obj.trustRegionRhoBound) * predictedDrop ...
                        && ~obj.candidateSolution.TrustRegionInterior
                    obj.trustRegionScaled = obj.SCALE_UP;
                else
                    obj.trustRegionScaled = obj.NO_SCALE;
                end
                
                if trueDrop > obj.trustRegionEta * predictedDrop
                    obj.newSolutionAccepted = true;
                else
                    obj.newSolutionAccepted = false;
                end
            else
                obj.trustRegionScaled = obj.NO_SCALE;
                obj.newSolutionAccepted = false;
            end
            
            obj.trustRegionRho = trueDrop / predictedDrop;
        end
        
        % UpdateBestSolution updates the best known solution, using the
        % true value of theta
        function UpdateBestSolution( obj )
            if obj.newSolutionAccepted
                if ~isempty( obj.bestSolution )
                    obj.secondBestSolution = obj.bestSolution.copy;
                end
                obj.bestSolution = obj.candidateSolution.copy;
                obj.CalculateProbability();
            end
        end
        
        % CalculateProbabilty calculates the worst case distribution for
        % the current best solution
        function CalculateProbability( obj )
            obj.pWorst = obj.bestSolution.Lambda*obj.numObsPerScen ...
                ./(obj.bestSolution.Mu-obj.bestSolution.SecondStageValues');
            obj.relativeLikelihood ...
                = exp(sum( obj.numObsPerScen.*( log(obj.pWorst) ...
                -log(obj.numObsPerScen/obj.numObsTotal) ) ));
        end
        
        % ResetSecondStageSolutions clears the second stage solution values
        % and dual solution information
        function ResetSecondStageSolutions( obj )
            obj.candidateSolution.Reset();
            obj.newSolutionAccepted = [];
            obj.zLowerUpdated = [];
            obj.trustRegionScaled = [];
        end
        
        % GetMasterc gets the cost vector for the master problem
        function cOut = GetMasterc( obj )
            cOut = [obj.lpModel.c, 0, 0, 0] * obj.objectiveScale;
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
            lOut(obj.LAMBDA) = obj.bestSolution.Lambda / 100;
            lOut(obj.MU) = -Inf;
            lOut(obj.THETA) = -Inf;
            
            obj.trustRegionLower = obj.GetDecisions(obj.bestSolution) ...
                - obj.trustRegionSize;
            obj.trustRegionLower(obj.THETA) = -Inf;
            
            lOut = max( lOut, obj.trustRegionLower );
        end
        
        % GetMasteru gets the upper bound for the master problem
        function uOut = GetMasteru( obj )
            uOut = [obj.lpModel.u; 0; 0; 0];
            uOut(obj.LAMBDA) = Inf;
            uOut(obj.MU) = Inf;
            uOut(obj.THETA) = Inf;
            
            obj.trustRegionUpper = obj.GetDecisions(obj.bestSolution) ...
                + obj.trustRegionSize;
            obj.trustRegionUpper(obj.THETA) = Inf;
            
            uOut = min( uOut, obj.trustRegionUpper );
        end
        
        % GetDecisions takes the first stage data from a Solution object
        % and formats it in a vector
        function vOut = GetDecisions( obj, solution, inType )
            if nargin < 3
                inType = 'true';
            end
            vOut = [solution.X; 0; 0; 0];
            vOut(obj.LAMBDA) = solution.Lambda;
            vOut(obj.MU) = solution.Mu;
            switch inType
                case 'master'
                    vOut(obj.THETA) = solution.ThetaMaster;
                case 'true'
                    vOut(obj.THETA) = solution.ThetaTrue;
                otherwise
                    error( 'Only accepts ''master and ''true''' )
            end
        end
        
    end

    %     Accessor methods
    methods (Access=public)
        % CandidateVector returns the current candidate vector
        function outCV = CandidateVector( obj )
            outCV = obj.GetDecisions( obj.candidateSolution, 'master' );
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
