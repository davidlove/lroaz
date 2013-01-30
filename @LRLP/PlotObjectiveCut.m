% Plot the cuts using the linear constraint A & b
function PlotObjectiveCut( obj, inVariableNumber, inCutNumber, inBounds )

A = obj.objectiveCutsMatrix;
b = obj.objectiveCutsRHS;
soln = obj.candidateSolution;
soln(obj.THETA) = 0;
cMaster = obj.GetMasterc();

slope = A( inCutNumber, inVariableNumber );
intercept = A( inCutNumber, [1:inVariableNumber-1, inVariableNumber+1:end-1]) ...
    * soln( [1:inVariableNumber-1, inVariableNumber+1:end-1] ) ...
    - b( inCutNumber );

x = inBounds([1,end]);
y = zeros(size(x));

for ii = 1:length(x)
    soln(inVariableNumber) = x(ii);
    y(ii) = intercept + slope*x(ii) + cMaster*soln;
end

if inCutNumber < obj.NumObjectiveCuts()
    color = 'g';
else
    color = 'm';
end

plot(x,y,color, 'LineWidth',1);

% rows = length(b);
% x = [-1 1];
% for ii=1:rows
%     slope = A(ii,1:end-3);
%     intercept = A(ii,end-2:end-1)*[lambda;mu] - b(ii);
%     y = [0 0];
%     y(1) = intercept + slope*x(1) + get_first_stage_obj(x(1),lambda,mu,c,N,Nbar);
%     y(2) = intercept + slope*x(2) + get_first_stage_obj(x(2),lambda,mu,c,N,Nbar);
%     plot(x,y,color, 'LineWidth',1);
% end
end