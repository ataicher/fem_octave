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
function [xPoints, yPoints, uExactVal] = evalExactSolution(appCtx)

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
    uExact = appCtx.field(f).uExact;
  
    uExactValTmp = zeros(numComp,numQuadPoints*numCells);
    for c=1:numCells
        realQuadPoints = projectQuadPoints(c,appCtx);    
        
        uExactValLocal = zeros(numComp,numQuadPoints);   
        for q=1:numQuadPoints
            for comp=1:numComp
                uExactValLocal(comp,q) = uExact{comp}(realQuadPoints(1,q),realQuadPoints(2,q));
            end
        end
        
        uExactValTmp(:,(c-1)*numQuadPoints+1:c*numQuadPoints) = uExactValLocal;
    end
    
    uExactVal.field{f} = uExactValTmp;
end
