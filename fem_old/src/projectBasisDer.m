% project the local degrees of freedom to function values at quadrature
% points
% INPUTS:
%   uLocal - local degrees of freedom
%   appCtx - problem parameters
% OUTPUT:
%   uLocalVal - uLocalVal(comp,f,:) is the function value of component comp of
%       field f
function realBasisDer = projectBasisDer(c,appCtx)

dim = appCtx.dim;
numFields = appCtx.numFields;
numQuadPoints = appCtx.quad.numQuadPoints;
J = appCtx.cellGeometry.J{c};
invJ = appCtx.cellGeometry.invJ{c};
detJ = appCtx.cellGeometry.detJ{c};
cells = appCtx.mesh.cells;


realBasisDer.field = cell(numFields,1);
for f=1:numFields
    numComp = appCtx.field(f).numComp;
    numBasisFuncs = appCtx.field(f).numBasisFuncs;
    basisInfo = appCtx.field(f).basisInfo;
    rBasisDer = appCtx.field(f).basisDer;
    
    for b=1:numBasisFuncs
        DOFType = basisInfo(b).DOFType;
        geoNum = basisInfo(b).geoNum;
        
        switch DOFType
                
            case 'normalFlux'
                
                lEPtr = geoNum;
                orient = sign(cells(lEPtr,c));
                for q=1:numQuadPoints
                    rBasisDer(b,:,:,q) = orient*invJ'*reshape(rBasisDer(b,:,:,q),dim,numComp)*J'/detJ;
                end
                
            case 'tangentFlux'
                
                lEPtr = geoNum;
                orient = sign(cells(lEPtr,c));
                for q=1:numQuadPoints
                    rBasisDer(b,:,:,q) = orient*invJ'*reshape(rBasisDer(b,:,:,q),dim,numComp)*invJ;
                end
                
            otherwise % 'nodal1' || 'nodal2' || 'nodal' || 'averageVal1' || 'averageVal2' || 'averageVal'     
                
                for q=1:numQuadPoints
                    rBasisDer(b,:,:,q) = invJ'*reshape(rBasisDer(b,:,:,q),dim,numComp);
                end
                
        end
    end
    
    realBasisDer.field{f} = rBasisDer;
end
