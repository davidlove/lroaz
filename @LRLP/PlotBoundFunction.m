function PlotBoundFunction( obj, inVariableNumber, inBounds )

A = obj.objectiveCutsMatrix;
b = obj.objectiveCutsRHS;
soln = obj.GetDecisions( obj.candidateSolution );
soln(obj.THETA) = 0;
cMaster = obj.GetMasterc();

x = linspace(inBounds(1),inBounds(end),101).';
theta = zeros(length(obj.THETA),length(x));

for thetaCount = 1:length(obj.THETA)
    rows = thetaCount:length(obj.THETA):(size(A,1)-1);
    
    slopes = A(rows,inVariableNumber);
    ints = A(rows,[1:inVariableNumber-1, inVariableNumber+1:end-1]) ...
        * soln([1:inVariableNumber-1, inVariableNumber+1:end-1]) ...
        - b(rows);
    
    sMat = repmat(slopes,1,length(x));
    iMat = repmat(ints,  1,length(x));
    xMat = repmat(x.',length(rows),1);
    
    assert( isequal( size(sMat), size(iMat) ) )
    
    thetas = sMat.*xMat + iMat;
    
    theta(thetaCount,:) = max(thetas,[],1);
end

soln( inVariableNumber ) = 0;
y = cMaster(obj.THETA)*theta + cMaster*soln + cMaster(inVariableNumber)*x.';

plot(x,y,'k--', 'LineWidth',2)