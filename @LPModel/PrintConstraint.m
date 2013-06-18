function PrintConstraint(obj, inConstraint, inPeriod)

if nargin < 3
    inPeriod = obj.timeLag + 1;
end

if ischar(inConstraint)
    constraintName = inConstraint;
    constraintRow = find(strcmp(constraintName,obj.Nr));
elseif isnumeric(inConstraint)
    constraintRow = inConstraint;
    constraintName = obj.Nr{constraintRow};
else
    error('inConstraint must be constraint name or index')
end

disp(['Constraint for ' constraintName ' in time period t = ' num2str(inPeriod) ':'])

constraint = '';

if inPeriod >= obj.timeLag + 1 % To include A_lag
    [~,c,v] = find(obj.A_lag(constraintRow,:));
    for ii=1:length(c)
        constraint = strcat( constraint, ...
            strsign(v(ii)), num2str(abs(v(ii))), ...
            '[', obj.variableNames(c(ii)), ']', ...
            '_[t-', num2str(obj.timeLag), ']' );
    end
end

if inPeriod >= 1 + 1 % To include A_st
    [~,c,v] = find(obj.A_st(constraintRow,:));
    for ii=1:length(c)
        constraint = strcat( constraint, ...
            strsign(v(ii)), num2str(abs(v(ii))), ...
            '[', obj.variableNames(c(ii)), ']', '_[t-1]' );
    end
end

[~,c,v] = find(obj.Abase(constraintRow,:));
for ii=1:length(c)
    constraint = strcat( constraint, ...
        strsign(v(ii)), num2str(abs(v(ii))), ...
        '[', obj.variableNames(c(ii)), ']', '_[t] ' );
end

btemp = reshape(obj.b,length(obj.Nr),obj.timePeriods);
constraint = strcat( constraint, ...
    ' = ', num2str(btemp(constraintRow,inPeriod)) );

disp(constraint{1})


function outSign = strsign(inNum)

switch sign(inNum)
    case {1,0}
        outSign = ' + ';
    case -1
        outSign = ' - ';
end