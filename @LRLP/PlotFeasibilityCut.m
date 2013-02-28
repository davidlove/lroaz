% Plot the cuts using the linear constraint A & b
function PlotFeasibilityCut( obj, inVariableNumber, inCutNumber, ylim )

A = obj.feasibilityCutsMatrix;
b = obj.feasibilityCutsRHS;
soln = obj.GetDecisions( obj.candidateSolution );
soln(obj.THETA) = 0;

slope = A( inCutNumber, inVariableNumber );
intercept = A( inCutNumber, [1:inVariableNumber-1, inVariableNumber+1:end-1]) ...
    * soln( [1:inVariableNumber-1, inVariableNumber+1:end-1] ) ...
    - b( inCutNumber );

if slope ~= 0
    y = ylim;
    x = -intercept/slope*ones(size(y));
    color = 'c';
    plot(x,y,color, 'LineWidth',1);
end

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