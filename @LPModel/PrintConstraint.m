function PrintConstraint(obj, inConstraint, inPeriod)

if nargin < 3
    inPeriod = obj.timeLag + 1;
end

if ischar(inConstraint)
    inConstraint = {inConstraint};
end

for ic = inConstraint
    if iscell(ic)
        constraintName = ic{1};
        constraintRow = find(strcmp(constraintName,obj.Nr));
        if isempty(constraintRow)
            error(['No constraint named "' constraintName '"'])
        end
    elseif isnumeric(ic)
        constraintRow = ic;
        constraintName = obj.Nr{constraintRow};
    else
        error('inConstraint must be constraint name or index')
    end
    
    disp(['Constraint for ' constraintName ' in time period t = ' num2str(inPeriod) ':'])
    
    constraint = strcat( WriteConstraint( obj, obj.A_lag, constraintRow, inPeriod, obj.timeLag ), ...
        WriteConstraint( obj, obj.A_st, constraintRow, inPeriod, 1 ), ...
        WriteConstraint( obj, obj.Abase, constraintRow, inPeriod, 0) );
    
    btemp = reshape(obj.b,length(obj.Nr),obj.timePeriods);
    constraint = sprintf( '%s%s%f', constraint, ...
        ' = ', num2str(btemp(constraintRow,inPeriod)) );
    
    disp(constraint)
    disp(' ')
end

function outString = WriteConstraint( obj, inMatrix, inRow, inPeriod, inDelay )
outString = '';
if inPeriod >= inDelay + 1 % Tests whether delay is great enough
    [~,c,v] = find(inMatrix(inRow,:));
    for ii=1:length(c)
        varName = obj.variableNames{c(ii)};
        varName(ismember(varName,' ')) = [];
        if inDelay == 0
            delayString = 't';
        else
            delayString = strcat( 't-', num2str(inDelay) );
        end
        outString = strcat( outString, sprintf( '%s%i[%s](%s)', ...
            strsign(v(ii)), abs(v(ii)), varName, delayString ) );
    end
end


function outSign = strsign(inNum)

switch sign(inNum)
    case {1,0}
        outSign = ' + ';
    case -1
        outSign = ' - ';
end