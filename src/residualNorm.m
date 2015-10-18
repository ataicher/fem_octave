% calculate the L2 error from the degrees of freedom using values at
% quadrature points
% INPUT:
%   u       - degrees of freedoms
%   appCtx  - problem specific information
% OUTPUT:
%   norm   - norm(field) is the L2 norm of the field
function norm = L2Norm(u, appCtx)

numCells = appCtx.mesh.numCells;
numQuadPoints = appCtx.quad.numQuadPoints;
quadWeights = appCtx.quad.quadWeights;
numFields = appCtx.numFields;

norm = zeros(numFields,1);
for c=1:numCells
    detJ = appCtx.cellGeometry.detJ{c};
    realBasis = projectBasis(c,appCtx);
    
    for f=1:numFields
        numComp = appCtx.field(f).numComp;
        cellDOF = appCtx.field(f).cellDOF(:,c);
        rBasis = realBasis.field{f};
        
        uLocal = u(cellDOF);
        for q=1:numQuadPoints
            
            for comp=1:numComp
                uVal = uLocal'*rBasis(:,comp,q);
                norm(f) = norm(f) + uVal^2*detJ*quadWeights(q);
            end
        end
    end
end
norm = sqrt(norm);

if sum(isnan(norm)) > 0
    error('residual norm cannot be calculated: NaN found')
end