% project a function to degrees of freedom of the AW basis functions
% INPUT:
%   fun - function to be projected. There are a total of numCompTotal
%       functions
%   appCtx - problem specific information
% OUTPUT:
%   u - degrees of freedom of projected function
function u = projectExactSolution(appCtx)

numFields = appCtx.numFields;
numCells = appCtx.mesh.numCells;
nodes = appCtx.mesh.nodes;
edges = appCtx.mesh.edges;
cells = appCtx.mesh.cells;
edgeNormals = appCtx.mesh.edgeNormals;
numQuadPoints = appCtx.quad.numQuadPoints;
numEdgeQuadPoints = appCtx.quad.numEdgeQuadPoints;
quadWeights = appCtx.quad.quadWeights;
edgeQuadWeights = appCtx.quad.edgeQuadWeights;
refNodes = appCtx.quad.refNodes;
refArea = appCtx.quad.refArea;
v_0Ref = appCtx.quad.v_0Ref;
globalSize = appCtx.globalSize;

u = zeros(globalSize,1);
for f=1:numFields
    numComp = appCtx.field(f).numComp;
    numBasisFuncs = appCtx.field(f).numBasisFuncs;
    basisInfo = appCtx.field(f).basisInfo;
    uExact = appCtx.field(f).uExact;
    
    for c = 1:numCells
        v_0 = appCtx.cellGeometry.v_0{c};
        J = appCtx.cellGeometry.J{c};
        detJ = appCtx.cellGeometry.detJ{c};
        cellDOF = appCtx.field(f).cellDOF(:,c);
        realQuadPoints = projectQuadPoints(c,appCtx);
        realArea = detJ*refArea;
        
        for b=1:numBasisFuncs
            DOFType = basisInfo(b).DOFType;
            geoNum = basisInfo(b).geoNum;
            
            switch DOFType
                case 'nodal1'
                    
                    n = v_0 + J*(refNodes(:,geoNum)-v_0Ref);
                    u(cellDOF(b)) = uExact{1}(n(1),n(2));
                    
                case 'nodal2'
                    
                    n = v_0 + J*(refNodes(:,geoNum)-v_0Ref);
                    u(cellDOF(b)) = uExact{2}(n(1),n(2));
                    
                case 'nodal'
                    
                    n = v_0 + J*(refNodes(:,geoNum)-v_0Ref);
                    u(cellDOF(b)) = uExact{1}(n(1),n(2));
                    
                case 'normalFlux'

                    realEdgeQuadPoints = projectEdgeQuadPoints(geoNum,c,appCtx);
                    eQuadWeights = edgeQuadWeights(:,geoNum);
                    ePtr = abs(cells(geoNum,c));
                    normal = edgeNormals(:,ePtr);
                    edge = nodes(:,edges(2,ePtr)) - nodes(:,edges(1,ePtr));
                    detJEdge = norm(J*edge)/norm(edge);
                    edgeInt = 0;
                    for q=1:numEdgeQuadPoints
                        for comp = 1:numComp
                            edgeInt = edgeInt + ...
                                detJEdge*normal(comp)*eQuadWeights(q)*uExact{comp}(realEdgeQuadPoints(1,q),realEdgeQuadPoints(2,q));  
                        end
                    end
                    u(cellDOF(b)) = edgeInt;
                    
                case 'tangentFlux'
                    
                    realEdgeQuadPoints = projectEdgeQuadPoints(geoNum,c,appCtx);
                    eQuadWeights = edgeQuadWeights(:,geoNum);
                    ePtr = abs(cells(geoNum,c));
                    tangent = [0 -1 ; 1 0]*edgeNormals(:,ePtr);
                    edge = nodes(:,edges(2,ePtr)) - nodes(:,edges(1,ePtr));
                    detJEdge = norm(J*edge)/norm(edge);
                    edgeInt = 0;
                    for q=1:numEdgeQuadPoints
                        for comp = 1:numComp
                            edgeInt = edgeInt + ...
                                detJEdge*tangent(comp)*eQuadWeights(q)*uExact{comp}(realEdgeQuadPoints(1,q),realEdgeQuadPoints(2,q));
                        end
                    end
                    u(cellDOF(b)) = edgeInt;
                    
                case 'averageVal1'
                    
                    int = 0;
                    for q=1:numQuadPoints
                        int = int + ...
                            (1/realArea)*detJ*quadWeights(q)*uExact{1}(realQuadPoints(1,q),realQuadPoints(2,q));
                    end
                    u(cellDOF(b)) = int;
                    
                case 'averageVal2'
                    
                    int = 0;
                    for q=1:numQuadPoints
                        int = int + ...
                            (1/realArea)*detJ*quadWeights(q)*uExact{2}(realQuadPoints(1,q),realQuadPoints(2,q));
                    end
                    u(cellDOF(b)) = int;
                    
                case 'averageVal'
                    
                    int = 0;
                    for q=1:numQuadPoints
                        int = int + ...
                            (1/realArea)*detJ*quadWeights(q)*uExact{1}(realQuadPoints(1,q),realQuadPoints(2,q));
                    end
                    u(cellDOF(b)) = int;
                    
            end
        end
    end
end

% test for NaN and replace with 0
u(isnan(u)) = 0;

