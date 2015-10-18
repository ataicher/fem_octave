% project the local degrees of freedom to function values at quadrature
% points
% INPUTS:
%   uLocal - local degrees of freedom
%   appCtx - problem parameters
% OUTPUT:
%   uLocalVal - uLocalVal(comp,f,:) is the function value of component comp of
%       field f
function realEdgeBasis = projectEdgeBasis(lEPtr,c,appCtx)

numFields = appCtx.numFields;
numEdgeQuadPoints = appCtx.quad.numEdgeQuadPoints;
J = appCtx.cellGeometry.J{c};
invJ = appCtx.cellGeometry.invJ{c};
detJ = appCtx.cellGeometry.detJ{c};
cells = appCtx.mesh.cells;

% obtain edge orientation
orient = sign(cells(lEPtr,c));

realEdgeBasis.field = cell(numFields,1);
for f=1:numFields
    numComp = appCtx.field(f).numComp;
    numBasisFuncs = appCtx.field(f).numBasisFuncs;
    basisInfo = appCtx.field(f).basisInfo;
    edgeBasis = reshape(appCtx.field(f).edgeBasis(:,:,:,lEPtr),numBasisFuncs,numComp,numEdgeQuadPoints);
    
    rEdgeBasis = zeros(numBasisFuncs,numComp,numEdgeQuadPoints);
    for b=1:numBasisFuncs
        DOFType = basisInfo(b).DOFType;
        
        switch DOFType
            case 'nodal1'
                
                rEdgeBasis(b,:,:) = edgeBasis(b,:,:);
                
            case 'nodal2'
                
                rEdgeBasis(b,:,:) = edgeBasis(b,:,:);
                
            case 'nodal'
                
                rEdgeBasis(b,1,:) = edgeBasis(b,1,:);
                
            case 'normalFlux'
                
                % realBasis = orient*J*basis
                for q=1:numEdgeQuadPoints
                    rEdgeBasis(b,:,q) = orient*edgeBasis(b,:,q)*J'/detJ;
                end
                
            case 'tangentFlux'
                
                % realBasis = orient*invJ'*basis
                for q=1:numEdgeQuadPoints
                    rEdgeBasis(b,:,q) = orient*edgeBasis(b,:,q)*invJ;
                end
                
            case 'averageVal1'
                
                rEdgeBasis(b,:,:) = edgeBasis(b,:,:);
                
            case 'averageVal2'
                
                rEdgeBasis(b,:,:) = edgeBasis(b,:,:);
                
            case 'averageVal'
                
                rEdgeBasis(b,1,:) = edgeBasis(b,1,:);
                         
        end
    end
    
    realEdgeBasis.field{f} = rEdgeBasis;
end

