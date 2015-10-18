% Project the degrees of freedom of a function to its values on quadrature
%   points
% INPUT:
%   u       - degrees of freedoms
%   appCtx  - problem specific information
% OUTPUT:
%   xPoints - x coordinate of points where the function is evaluated
%   yPoints - y coordinate of points where the function is evaluated
%   funcVal - funcVal(comp,f,:) is the function value of component comp of
%       field f
function [xPoints, yPoints, funDerVal] = projectDOFDeriv(u, appCtx)
dim = appCtx.dim;
quad = appCtx.quad;
mesh = appCtx.mesh;
variables = appCtx.variables;
LG = appCtx.LG;

numQuadPoints = quad.numQuadPoints;
numCells = mesh.numCells;
numFields = appCtx.variables.numFields;

xPoints = zeros(numQuadPoints*numCells,1); yPoints = zeros(numQuadPoints*numCells,1);
for c=1:numCells
    realQuadPoints = projectQuadPoints(c,appCtx);
    xPoints((c-1)*numQuadPoints+1:c*numQuadPoints) = realQuadPoints(1,:)';
    yPoints((c-1)*numQuadPoints+1:c*numQuadPoints) = realQuadPoints(2,:)';
end

for f=1:numFields
    numComp = variables.field{f}.numComp;
    cellDOF = LG.field{f}.cellDOF;
    
    funDerVal.field{f} = zeros(dim,numComp,numQuadPoints*numCells);
    for c=1:numCells
        uLocal = u(cellDOF(:,c));
        realBasisDer = projectBasisDer(c,appCtx);
        rBasisDer = realBasisDer.field{f};

        funDerValLocal = zeros(dim,numComp,numQuadPoints);
        for q=1:numQuadPoints
            for comp=1:numComp
                for d=1:dim
                    funDerValLocal(d,comp,q) = uLocal'*rBasisDer(:,d,comp,q);
                end
            end
        end
        
        funDerVal.field{f}(:,:,(c-1)*numQuadPoints+1:c*numQuadPoints) = funDerValLocal;
    end
end
