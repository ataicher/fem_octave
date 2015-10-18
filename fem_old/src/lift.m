%% LIFT
% generate an initial guess u that matches the given essential conditions.
% Project essential condtion function to DOF's on boundary
%
function u = lift(appCtx)

numFields = appCtx.numFields;
nodes = appCtx.mesh.nodes;
edges = appCtx.mesh.edges;
cells = appCtx.mesh.cells;
numBndryEdges = appCtx.mesh.numBndryEdges;
bndryEdges = appCtx.mesh.bndryEdges;
edgeCells = appCtx.mesh.edgeCells;
edgeNormals = appCtx.mesh.edgeNormals;
numLocalEdges = appCtx.mesh.numLocalEdges;
numEdgeQuadPoints = appCtx.quad.numEdgeQuadPoints;
edgeQuadWeights = appCtx.quad.edgeQuadWeights;
globalSize = appCtx.globalSize;

u = zeros(globalSize,1);
for f=1:numFields
    
    % check if boundary condition exists
    if appCtx.field(f).BNDRYCONDTION
        numComp = appCtx.field(f).numComp;
        numBasisFuncs = appCtx.field(f).numBasisFuncs;
        basisInfo = appCtx.field(f).basisInfo;
        edgeBC = appCtx.field(f).edgeBC;
           
        for bEPtr=1:numBndryEdges
            bndry = appCtx.field(f).bndry(edgeBC(bEPtr));
            alpha = bndry.alpha;
            beta = bndry.beta;
            eta = bndry.eta;
            
            % determine the cell, edge orientation, and local edge
            ePtr = bndryEdges(bEPtr);
            c = edgeCells(1,ePtr);
            lEPtr = find(abs(cells(:,c))==ePtr);
            orient = sign(cells(lEPtr,c));
            cellDOF = appCtx.field(f).cellDOF(:,c);
            
            % get normal and tangent direction on edge
            normal = orient*edgeNormals(:,ePtr);
            tangent = [0 -1 ; 1 0]*normal;
            
            % get global nodes and corresponding local nodes
            n1Ptr = edges(1,ePtr);
            n2Ptr = edges(2,ePtr);
            n1 = nodes(:,n1Ptr);
            n2 = nodes(:,n2Ptr);
            if orient == 1
                lN1Ptr = lEPtr;
                lN2Ptr = mod(lEPtr,numLocalEdges) + 1;
            else
                lN1Ptr = mod(lEPtr,numLocalEdges) + 1;
                lN2Ptr = lEPtr;
            end
            
            % get edge quadrature and determinant to evaluate edge integral
            realQuadPoints = projectEdgeQuadPoints(lEPtr,c,appCtx);
            realEdgeQuadWeights = edgeQuadWeights(:,lEPtr);
            J = appCtx.cellGeometry.J{c};
            e = n2-n1;
            detJEdge = norm(J*e)/norm(e);
            
            % scalar field
            if numComp == 1
                
                if alpha == 0
                    for b = 1:numBasisFuncs
                        DOFType = basisInfo(b).DOFType;
                        geoNum = basisInfo(b).geoNum;
                        
                        switch DOFType
                            case 'nodal'
                                if geoNum == lN1Ptr
                                    u(cellDOF(b)) = eta{1}(n1(1),n1(2))/beta;
                                elseif geoNum == lN2Ptr
                                    u(cellDOF(b)) = eta{1}(n2(1),n2(2))/beta;
                                end
                        end
                    end
                end
                
            % vector field
            else
                for b = 1:numBasisFuncs
                    DOFType = basisInfo(b).DOFType;
                    geoNum = basisInfo(b).geoNum;
                    
                    switch DOFType
                        case 'nodal1'
                            
                            if geoNum == lN1Ptr
                                u(cellDOF(b)) = (alpha(1)==0)*eta{1}(n1(1),n1(2))*normal(1)/(beta(1) + (beta(1)==0)) + ...;
                                    + (alpha(2)==0)*eta{2}(n1(1),n1(2))*tangent(1)/(beta(2) + (beta(2)==0));
                            elseif geoNum ==lN2Ptr
                                u(cellDOF(b)) = (alpha(1)==0)*eta{1}(n2(1),n2(2))*normal(1)/(beta(1) + (beta(1)==0)) + ...;
                                    + (alpha(2)==0)*eta{2}(n2(1),n2(2))*tangent(1)/(beta(2) + (beta(2)==0));
                            end
                            
                        case 'nodal2'
                            
                            if geoNum == lN1Ptr
                                u(cellDOF(b)) = (alpha(1)==0)*eta{1}(n1(1),n1(2))*normal(2)/(beta(1) + (beta(1)==0)) + ...;
                                    + (alpha(2)==0)*eta{2}(n1(1),n1(2))*tangent(2)/(beta(2) + (beta(2)==0));
                                
                            elseif geoNum == lN2Ptr
                                u(cellDOF(b)) = (alpha(1)==0)*eta{1}(n2(1),n2(2))*normal(2)/(beta(1) + (beta(1)==0)) + ...;
                                    + (alpha(2)==0)*eta{2}(n2(1),n2(2))*tangent(2)/(beta(2) + (beta(2)==0));
                            end
                            
                        case 'normalFlux'
                            
                            if geoNum == lEPtr
                                if alpha(1) == 0
                                    for q=1:numEdgeQuadPoints
                                        u(cellDOF(b)) = u(cellDOF(b)) + ...
                                            orient*detJEdge*realEdgeQuadWeights(q)*eta{1}(realQuadPoints(1,q),realQuadPoints(2,q))/beta(1);
                                    end
                                end
                            end
                            
                        case 'tangentFlux'
                            
                            if geoNum == lEPtr
                                if alpha(2) == 0
                                    for q=1:numEdgeQuadPoints
                                        u(cellDOF(b)) = u(cellDOF(b)) + ...
                                            orient*detJEdge*realEdgeQuadWeights(q)*eta{2}(realQuadPoints(1,q),realQuadPoints(2,q))/beta(2);
                                    end
                                end
                            end
                    end
                end
            end
        end
    end
    
    % eliminate constant nullspace
    if appCtx.field(f).NULLSPACE
        nullDOF = appCtx.field(f).NULLSPACE;
        if appCtx.EXISTEXACTSOL
            % project exact solution to degrees of freedom
            uExactSol = projectExactSolution(appCtx);
            u(nullDOF) = uExactSol(nullDOF);
        else
            u(nullDOF) = 0;
        end
    end
    
end
