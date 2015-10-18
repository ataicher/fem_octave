% calculate the L2 error from the degrees of freedom using values at
% quadrature points
% INPUT:
%   u       - degrees of freedoms
%   appCtx  - problem specific information
% OUTPUT:
%   norm   - norm(field) is the L2 norm of the field
function norm = L2NormPhi(u, appCtx)

numCells = appCtx.mesh.numCells;
numQuadPoints = appCtx.quad.numQuadPoints;
quadWeights = appCtx.quad.quadWeights;
numFields = appCtx.numFields;
phi = appCtx.phi;

for f = 1:numFields
    name = appCtx.field(f).name;
    if strcmp(name,'fluid pressure potential')
        fFluidPres = f;
    end
end

norm = 0;
for c=1:numCells
    detJ = appCtx.cellGeometry.detJ{c};
    realBasis = projectBasis(c,appCtx);
    realQuadPoints = projectQuadPoints(c,appCtx);
    
        numComp = appCtx.field(fFluidPres).numComp;
        cellDOF = appCtx.field(fFluidPres).cellDOF(:,c);
        rBasis = realBasis.field{fFluidPres};
        name = appCtx.field(fFluidPres).name;
        
        uLocal = u(cellDOF);
            
            for q=1:numQuadPoints
                phiVal = phi(realQuadPoints(1,q),realQuadPoints(2,q));
                for comp=1:numComp
                    
                    uVal = sqrt(phiVal)*uLocal'*rBasis(:,comp,q);
                    norm(f) = norm(f) + uVal^2*detJ*quadWeights(q);
                end
            end
            
        else
            
            for q=1:numQuadPoints
                
                for comp=1:numComp
                    uVal = uLocal'*rBasis(:,comp,q);
                    norm(f) = norm(f) + uVal^2*detJ*quadWeights(q);
                end
            end
        end
    
end
norm = sqrt(norm);