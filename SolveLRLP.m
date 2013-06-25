function [lrlp, outTotalCuts, outTotalProbs] = ...
    SolveLRLP( simpleLP, phi, obs, rho, cutType, isGraphical )

if nargin < 6
    isGraphical = false;
    if nargin < 5
        cutType = 'multi';
    end
end

if ~ismember( cutType, {'single','multi'} )
        error([cutType ' is not a valid cut type'])
end

opt = 'cplexlp';
if exist(opt,'file')
    % Keep default optimizer
else
    disp(' ')
    disp(['Optimizer ' opt ' not found.'])
    opt = 'linprog';
    disp(['Falling back to ' opt '.'])
    disp(' ')
end

lrlp = LRLP( simpleLP, phi, obs, rho, opt, cutType );

totalProblemsSolved = 1;
totalCutsMade = 1;
while lrlp.currentObjectiveTolerance > lrlp.objectiveTolerance
    totalProblemsSolved = totalProblemsSolved + 1;
    
    exitFlag = lrlp.SolveMasterProblem;
    if exitFlag ~= 1 || (exist('cS','var') && isequal(cS,lrlp.CandidateVector))
        if cS == lrlp.CandidateVector
            exitFlag = -100;
        end
        disp(['exitFlag = ' num2str(exitFlag)])
        switch exitFlag
            case 0
                lrlp.DoubleIterations;
            case {-2,-3,-4,-5}
                % Do nothing extra
            case -50
                % The optimizer failed to find a solution better than
                % lrlp.bestSolution
                % This has been observed with linprog and cplexlp
            case -100
                % The optimizer returned the same solution as it found the
                % previous time around
                % This has been observed with cplexlp
            otherwise
                error( ['Unknown error code: ' num2str(exitFlag)] )
        end
        lrlp.DeleteOldestCut;
        if lrlp.NumFeasibilityCuts > 1
            lrlp.DeleteOldestFeasibilityCut;
        end
        disp([num2str(lrlp.NumObjectiveCuts) ' Objective Cuts Remaining, ' ...
            num2str(lrlp.NumFeasibilityCuts) ' Feasibility Cuts Remaining.'])
        continue
    end
    cS = lrlp.CandidateVector;
    
    lrlp.SolveSubProblems;
    
    totalCutsMade = totalCutsMade + 1;
    lrlp.GenerateCuts;
    
    % Note, putting plot after updating trust region or solution will
    % result in it plotting the wrong trust region
    if isGraphical
        lrlp.Plot;
    end
    
    lrlp.UpdateTrustRegionSize;
    
    lrlp.UpdateSolutions;
    
    lrlp.UpdateTolerances;
    
    lrlp.WriteProgress;
    
    disp(['Total cuts made: ' num2str(totalCutsMade)])
    disp(['Total problems solved: ' num2str(totalProblemsSolved)])
    
    if isGraphical
        pause(.5)
    end
end

outTotalCuts = totalCutsMade;
outTotalProbs = totalProblemsSolved;

end
