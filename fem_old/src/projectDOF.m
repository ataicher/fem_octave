% Project the degrees of freedom of a function to its values on quadrature
%   points
% INPUT:
%   u       - degrees of freedoms
%   appCtx  - problem specific information
% OUTPUT:
%   xPoints - x coordinate of points where the function is evaluated
%   yPoints - y coordinate of points where the function is evaluated
%   funVal  - funVal(comp,f,:) is the function value of component comp of
%       field f
function [xPoints, yPoints, funVal] = projectDOF(u, appCtx)

numFields = appCtx.numFields;
numQuadPoints = appCtx.quad.numQuadPoints;
numCells = appCtx.mesh.numCells;

xPoints = zeros(numQuadPoints*numCells,1); yPoints = zeros(numQuadPoints*numCells,1);
for c=1:numCells
    realQuadPoints = projectQuadPoints(c,appCtx);
    xPoints((c-1)*numQuadPoints+1:c*numQuadPoints) = realQuadPoints(1,:)';
    yPoints((c-1)*numQuadPoints+1:c*numQuadPoints) = realQuadPoints(2,:)';
end
          
for f=1:numFields
    numComp = appCtx.field(f).numComp;
    cellDOF = appCtx.field(f).cellDOF;
    
    funVal.field{f} = zeros(numComp,numQuadPoints*numCells);
    for c=1:numCells
        uLocal = u(cellDOF(:,c));
        realBasis = projectBasis(c,appCtx);
        rBasis = realBasis.field{f};
            
        funValLocal = zeros(numComp,numQuadPoints);   
        for q=1:numQuadPoints
            for comp=1:numComp
                funValLocal(comp,q) = uLocal'*rBasis(:,comp,q);
            end
        end
        
        funVal.field{f}(:,(c-1)*numQuadPoints+1:c*numQuadPoints) = funValLocal;
    end
end
