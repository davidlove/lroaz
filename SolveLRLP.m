function [lrlp, outTotalCuts, outTotalProbs] = SolveLRLP( varargin )
% SolveLRLP implements the Bender's Decomposition method to solve the
% Phi-divergence problem
%
% Keyword arguments (*required arguments):
%    cuttype: 'single' cut or 'multi' cut
%    isgraphical: true or false for desired graphical optimization (very slow)
%    lp*: LPModel object defining the SLP-2
%    obs*: vector of observations numbers of the scenarios
%    phi*: PhiDivergence object or string defining the phi-divergence
%    rho*: Value of rho to construct the confidence region

requiredArgs = {'lp', 'obs', 'phi', 'rho'};

if mod(length(varargin),2) == 1
    error('Arguments must be key, value pairs')
end

for vv = 1:2:length(varargin)
    key = varargin{vv};
    value = varargin{vv+1};
    switch key
        case 'cuttype'
            cutType = value;
        case 'isgraphical'
            isGraphical = value;
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
        case 'rho'
            rho = value;
        otherwise
            error(['Unknown variable ', key])
    end
end

if ~exist('cutType', 'var')
    cutType = 'multi';
end
if ~exist('isGraphical', 'var')
    isGraphical = false;
end

if ~ismember( cutType, {'single','multi'} )
        error([cutType ' is not a valid cut type'])
end

for aa = requiredArgs
    if ~exist(aa{1}, 'var')
        error([aa{1}, ' is required but not defined'])
    end
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

lrlp = LRLP( lp, phi, obs, rho, opt, cutType );

totalProblemsSolved = 1;
totalCutsMade = 1;
while lrlp.currentObjectiveTolerance > lrlp.objectiveTolerance
    totalProblemsSolved = totalProblemsSolved + 1;
    
    exitFlag = lrlp.SolveMasterProblem;
    if (exist('cS','var') && isequal(cS,lrlp.CandidateVector))
        disp(' ')
        disp('Repeat Solution')
        exitFlag = -100;
    end
    if exitFlag ~= 1 && exitFlag ~= -100
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
%         if lrlp.NumFeasibilityCuts > 1
            lrlp.DeleteOldestFeasibilityCut;
%         end
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
    
    if exitFlag == -100
        lrlp.ForceAcceptSolution();
    end
    
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
