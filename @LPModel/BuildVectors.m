function [b1,UB1,LB1,Cost1,b2,UB2,LB2,Cost2,costOut] = BuildVectors(obj, cellInputFile)

% BuildVectors builds all vectors (right hand side b, upper bound UB, lower
% bound LB, cost vector Cost) for both the first and second stages.  First
% stage vectors names have a 1 at the end, second stage vector names have a
% 2 at the end.

numFiles = length(cellInputFile);

% Read input files into three-dimensional arrays
costVector = ReadInput(obj, cellInputFile, 'Costs');
lb = ReadInput(obj, cellInputFile, 'Lower Bounds');
ub = ReadInput(obj, cellInputFile, 'Upper Bounds');
RHS = ReadInput(obj, cellInputFile, 'b_vec');

f = find(ub==-1);
ub(f) = Inf(size(f));

if obj.timeLag ~= 1
    assert(obj.timeLag == 12);
    
    numVariablesPerPeriod = size(costVector,1);
    numConstraintsPerPeriod = size(RHS,1);
    
    assert( size(lb,1) == numVariablesPerPeriod );
    assert( size(ub,1) == numVariablesPerPeriod );
    
    costExp = zeros( numVariablesPerPeriod, obj.timePeriods, numFiles );
    rhsExp = zeros( numConstraintsPerPeriod, obj.timePeriods, numFiles );
    lbExp = zeros( numVariablesPerPeriod, obj.timePeriods, numFiles );
    ubExp = zeros( numVariablesPerPeriod, obj.timePeriods, numFiles );
    
    for ff = 1:numFiles
        sf = xlsread( strrep(cellInputFile{ff},'Inputs','sf'),'sf' );
        
        for jj = 1:obj.numYears
            for ii = 1:obj.timeLag
                kk = (jj-1)*obj.timeLag + ii;
                costExp(:,kk,ff) = costVector(:,jj,ff);
                ubExp(:,kk,ff) = ub(:,jj,ff) / obj.timeLag;
                lbExp(:,kk,ff) = lb(:,jj,ff) / obj.timeLag;
                
                rhsExp(1:15,kk,ff) = RHS(1:15,jj,ff) / obj.timeLag;
                rhsExp(16:end,kk,ff) = RHS(16:end,jj,ff) .* sf(ii,1);
                rhsExp(64:100,kk,ff) = RHS(64:100,jj,ff) / obj.timeLag;
            end
        end
        costVector = costExp;
        RHS = rhsExp;
        lb = lbExp;
        ub = ubExp;
    end
end

if obj.timePeriods > min([size(costVector,2),size(ub,2),size(lb,2),size(RHS,2)])
    error('Request for more periods than data provides');
end

% Remove unnecessary time periods
costVector = costVector(:,1:obj.timePeriods,:);
ub = ub(:,1:obj.timePeriods,:);
lb = lb(:,1:obj.timePeriods,:);
RHS = RHS(:,1:obj.timePeriods,:);

% Decompose into two stages
[b1,b2] = DecomposeStages(obj, RHS, true);
[UB1, UB2] = DecomposeStages(obj, ub, true);
[LB1, LB2] = DecomposeStages(obj, lb, false);
[Cost1, Cost2] = DecomposeStages(obj, costVector, false);

costOut = costVector(:,:,1);

function outMatrix = ReadInput( obj, cellInputFile, sheet )
% ReadInput reads Inputs.xlsx into three-dimensional arrays

numFiles = length(cellInputFile);

switch sheet
    case {'Costs', 'Upper Bounds', 'Lower Bounds'}
        prediction = obj.variableNames;
    case 'b_vec'
        prediction = obj.Nr;
    otherwise
        error(['Unassigned sheet name ', sheet])
end

for ii=1:numFiles
    [inMat,txt] = xlsread(cellInputFile{ii},sheet);
    if ~isequal( prediction, ...
            txt(2:2+length(prediction)-1,2) )
        num = length(prediction);
        errors = false(1,num);
        for nn = 1:num
            errors(nn) = ~isequal( prediction{nn}, txt{1+nn,2} );
        end
        errnums = find(errors)';
        errlist = [{'Excel row', 'Predicted name', 'Excel name'}; ...
            num2cell(errnums+1), prediction(errors), txt(1+errnums,2)];
        fprintf('\nERROR PRINTOUT')
        disp(errlist)
        error(['Predicted variable names in ''', sheet, ''' do not match those in file', cellInputFile{ii}])
    end
    if ii == 1
        outMatrix = zeros(size(inMat,1)-1, size(inMat,2), numFiles);
    end
    outMatrix(:,:,ii) = inMat(2:end,:);
end

function [stageOne, stageTwo] = DecomposeStages(obj, matrix, adjust)
% DecomposeStages splits the matrix into vectors for each stage

numFiles = size(matrix,3);
for ii=1:numFiles
    stage1Temp = matrix(:,1:obj.firstStagePeriods,ii);
    if ii == 1
        stageOne = zeros(numel(stage1Temp), numFiles);
    end
    stageOne(:,ii) = stage1Temp(:);
    stage2Temp = matrix(:,obj.firstStagePeriods+1:end,ii);
    if ii > 1 && adjust
        for jj=1:size(stage2Temp,2)
            stage2Temp(:,jj) = stage2Temp(:,jj) ...
                - matrix(:,obj.firstStagePeriods,ii) + matrix(:,obj.firstStagePeriods,1);
        end
    end
    if ii == 1
        stageTwo = zeros(numel(stage2Temp), numFiles);
    end
    stageTwo(:,ii) = stage2Temp(:);
end

stageOne = stageOne(:,1);
