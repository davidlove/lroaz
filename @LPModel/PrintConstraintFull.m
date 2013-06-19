function PrintConstraintFull(obj, inConstraint, inPeriod)

% PrintConstraintFull prints the requested constraint(s) in a human 
% readable format, and expands the loss variables automatically.
%
% PrintConstraint can be called in several ways:
% 1. obj.PrintConstraintFull( row ) prints the constraint by the specified
%    row number
% 2. obj.PrintConstraintFull( vector ) prints the constraints specified by
% each row number in the vector.
% 3. obj.PrintConstraintFull( string ) prints the constraint for the node
%    named in the input string
% 4. obj.PrintConstraintFull( cellArray ) prints the constraints for each
%    node named in the cell array of strings.
%
% An optional numeric argument specifies the time period to print the
% constraint.  If omitted, the earliest time period for which all lag and
% storage variables appear is chosen.

if nargin < 3
    inPeriod = obj.timeLag + 1;
end

if ischar(inConstraint)
    inConstraint = {inConstraint};
end

rhs = reshape(obj.b,length(obj.Nr),obj.timePeriods);

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
    
    disp(' ')
    disp(['Constraint for ' constraintName ' in time period t = ' num2str(inPeriod) ':'])
    
    constraint = strcat( WriteConstraint( obj, obj.A_lag, constraintRow, inPeriod, obj.timeLag ), ...
        WriteConstraint( obj, obj.A_st, constraintRow, inPeriod, 1 ), ...
        WriteConstraint( obj, obj.Abase, constraintRow, inPeriod, 0) );
    
    constraint = sprintf( '%s = %f', constraint, ...
        rhs(constraintRow,inPeriod) );
    
    disp(constraint)
end
disp(' ')

function outString = WriteConstraint( obj, inMatrix, inRow, inPeriod, inDelay )
numNodes = length(find(not(cellfun('isempty', strfind(obj.Nr, 'oss')))));
numNodes = size(obj.Nr,1) - numNodes;
outString = '';
if inPeriod >= inDelay + 1 % Tests whether delay is great enough
    [~,c,v] = find(inMatrix(inRow,:));
    for ii=1:length(c)
        varName = obj.variableNames{c(ii)};
        varName(ismember(varName,' ')) = [];
        if inRow+numNodes <= size(inMatrix,1)
            loss = full(inMatrix(inRow+numNodes,c(ii)));
        else
            loss = 0;
        end
        if inDelay == 0
            delayString = 't';
        else
            delayString = strcat( 't-', num2str(inDelay) );
        end
        if loss ~= 0
            if v(ii) + loss ~= 0
                outString = strcat( outString, sprintf( ' %s (%0.2f %s %0.2f)[%s](%s)', ...
                    strsign(v(ii)), abs(v(ii)), ...
                    strsign(v(ii)*loss), abs(loss), ...
                    varName, delayString ) );
            end
        else
            outString = strcat( outString, sprintf( ' %s %0.2f[%s](%s)', ...
                strsign(v(ii)), abs(v(ii)), varName, delayString ) );
        end
    end
end


function outSign = strsign(inNum)

switch sign(inNum)
    case {1,0}
        outSign = '+';
    case -1
        outSign = '-';
end