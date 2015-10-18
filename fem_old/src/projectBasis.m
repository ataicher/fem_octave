% project the local degrees of freedom to function values at quadrature
% points
% INPUTS:
%   uLocal - local degrees of freedom
%   appCtx - problem parameters
% OUTPUT:
%   uLocalVal - uLocalVal(comp,f,:) is the function value of component comp of
%       field f
function realBasis = projectBasis(c,appCtx)

numFields = appCtx.numFields;
numQuadPoints = appCtx.quad.numQuadPoints;
J = appCtx.cellGeometry.J{c};
invJ = appCtx.cellGeometry.invJ{c};
detJ = appCtx.cellGeometry.detJ{c};
cells = appCtx.mesh.cells;

realBasis.field = cell(numFields,1);
for f=1:numFields
    numBasisFuncs = appCtx.field(f).numBasisFuncs;
    basisInfo = appCtx.field(f).basisInfo;
    rBasis = appCtx.field(f).basis;
    
    for b=1:numBasisFuncs
        DOFType = basisInfo(b).DOFType;
        geoNum = basisInfo(b).geoNum;
        
        switch DOFType
               
            case 'normalFlux'
                
                % rBasis = (orient/detJ)*J*rBasis;
                lEPtr = geoNum;
                orient = sign(cells(lEPtr,c));
                for q=1:numQuadPoints
                    rBasis(b,:,q) = orient*rBasis(b,:,q)*J'/detJ;
                end
                
            case 'tangentFlux'
                
                % rBasis = orient*invJ'*rBasis;
                lEPtr = geoNum;
                orient = sign(cells(lEPtr,c));
                for q=1:numQuadPoints
                    rBasis(b,:,q) = orient*rBasis(b,:,q)*invJ;
                end
                
            otherwise % 'nodal1' || 'nodal2' || 'nodal' || 'averageVal1' || 'averageVal2' || 'averageVal' 
                
        end
    end
    
    realBasis.field{f} = rBasis;
end
