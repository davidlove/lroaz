function PrintConstraint(obj, inConstraint, varargin)

% PrintConstraint prints the requested constraint(s) in a human readable
% format.
%
% PrintConstraint can be called in several ways:
% 1. obj.PrintConstraint( row ) prints the constraint by the specified row
%    number
% 2. obj.PrintConstraint( vector ) prints the constraints specified by each
%    row number in the vector.
% 3. obj.PrintConstraint( string ) prints the constraint for the node named
%    in the input string
% 4. obj.PrintConstraint( cellArray ) prints the constraints for each node
%    named in the cell array of strings.
%
% An optional numeric argument specifies the time period to print the
% constraint.  If omitted, the earliest time period for which all lag and
% storage variables appear is chosen.
%
% An optional string argument 'loss' may be given, which will expand the
% loss variable and represent the losses as part of the constant on
% incoming arcs.

inPeriod = obj.timeLag + 1;
withLoss = false;

for ii = 1:length(varargin)
    if isnumeric(varargin{ii})
        inPeriod = varargin{ii};
    elseif ischar(varargin{ii}) && strcmpi(varargin{ii}, 'loss')
        withLoss = true;
    end
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
    
    constraint = strcat( WriteConstraint( obj, obj.A_lag, constraintRow, inPeriod, obj.timeLag, withLoss ), ...
        WriteConstraint( obj, obj.A_st, constraintRow, inPeriod, 1, withLoss ), ...
        WriteConstraint( obj, obj.Abase, constraintRow, inPeriod, 0, withLoss ) );
    
    constraint = sprintf( '%s = %f', constraint, ...
        rhs(constraintRow,inPeriod) );
    
    disp(constraint)
end
disp(' ')

function outString = WriteConstraint( obj, inMatrix, inRow, inPeriod, inDelay, inLoss )
numNodes = length(find(not(cellfun('isempty', strfind(obj.Nr, 'oss')))));
numNodes = size(obj.Nr,1) - numNodes;
outString = '';
if inPeriod >= inDelay + 1 % Tests whether delay is great enough
    [~,c,v] = find(inMatrix(inRow,:));
    for ii=1:length(c)
        varName = obj.variableNames{c(ii)};
        varName(ismember(varName,' ')) = [];
        if inLoss && inRow+numNodes <= size(inMatrix,1)
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