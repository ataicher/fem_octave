function [xPoints, yPoints, streamFun] = computeStreamFunction(u,appCtx)

dim = appCtx.dim;
numFields = appCtx.numFields;
numNodes = appCtx.mesh.numNodes;
nodes = appCtx.mesh.nodes;
numEdges = appCtx.mesh.numEdges;
edgeCells = appCtx.mesh.edgeCells;
edgeNormals = appCtx.mesh.edgeNormals;
edges = appCtx.mesh.edges;
cells = appCtx.mesh.cells;
edgeQuadWeights = appCtx.quad.edgeQuadWeights;
numEdgeQuadPoints = appCtx.quad.numEdgeQuadPoints;

xPoints = nodes(1,:);
yPoints = nodes(2,:);

edgeConnectivity = zeros(numEdges,numNodes);
n1Ptr = edges(1,1);
edgeConnectivity(1,n1Ptr) = 1;
for ePtr = 2:numEdges
    n1Ptr = edges(1,ePtr); n2Ptr = edges(2,ePtr);
    edgeConnectivity(ePtr,n1Ptr) = -1; edgeConnectivity(ePtr,n2Ptr) = 1;
end

streamFun = cell(numFields,1);
for f=1:numFields
    numComp = appCtx.field(f).numComp;
    
    if numComp == dim
        edgeFlux = zeros(numEdges,1);
        for ePtr = 1:numEdges
            c = edgeCells(1,ePtr);
            lEPtr = find(abs(cells(:,c))==ePtr);
            normal = edgeNormals(:,ePtr);
            uLocal = u(appCtx.field(f).cellDOF(:,c));
            realEdgeBasis = projectEdgeBasis(lEPtr,c,appCtx);
            rEdgeBasis = realEdgeBasis.field{f};
            J = appCtx.cellGeometry.J{c};
            edge = nodes(:,edges(2,ePtr)) - nodes(:,edges(1,ePtr));
            detJEdge = norm(J*edge)/norm(edge);

            for q=1:numEdgeQuadPoints

                normalComp = uLocal'*rEdgeBasis(:,:,q)*normal;               
                edgeFlux(ePtr) = edgeFlux(ePtr) + detJEdge*edgeQuadWeights(q)*normalComp;
            end
            
        end

        streamFun{f} = edgeConnectivity\edgeFlux;
    else
        streamFun{f} = 'NULL';
    end
end

