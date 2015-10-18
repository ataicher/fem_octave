%% COMPUTE CELL GEOMETRY---------------------------------------------------
% compute the transformation information for every cell in the mesh
%
% INPUT:
%   appCtx: struct containing problem information
%
% OUTPUT: cellGeometry struct containing
%   v_0:  bottom left vertex coordinates of element
%   J:    Jacobian of the linear transformation from reference element to 
%         real element
%   invJ: inverse Jacobian of the linear transformation from reference 
%         element to real element
%   detJ: determinant of the Jacobian of the linear transformation from 
%         reference element to real element
%
function cellGeometry = assembleCellGeometry(appCtx)

dim = appCtx.dim;
cells = appCtx.mesh.cells;
numCells = appCtx.mesh.numCells;
numLocalEdges = appCtx.mesh.numLocalEdges;
edges = appCtx.mesh.edges;
nodes = appCtx.mesh.nodes;
refNodes = appCtx.quad.refNodes;
v_0Ref = appCtx.quad.v_0Ref;


v_0 = cell(numCells,1);
J = cell(numCells,1);
invJ = cell(numCells,1);
detJ = cell(numCells,1);
for c=1:numCells
    
    % obtain nodes of real element
    realNodes = getRealNodes(c, dim, numLocalEdges, nodes, edges, cells);
    
    % choose 2 corresponding edges in the reference and real element
    e1Ref = refNodes(:,2)-v_0Ref;
    e2Ref = refNodes(:,3)-v_0Ref;
    ref = [e1Ref e2Ref];
    e1 = realNodes(:,2)-realNodes(:,1);
    e2 = realNodes(:,3)-realNodes(:,1);
    real = [e1 e2];

    v_0{c} = realNodes(:,1);
    J{c} = real*inv(ref);
    invJ{c} = inv(J{c});
    detJ{c} = det(J{c});    
end

    % collect into cellGeometry struct
    cellGeometry.v_0 = v_0;
    cellGeometry.J = J;
    cellGeometry.invJ = invJ;
    cellGeometry.detJ = detJ;    

end

%% GETREALNODES
% get the nodes of the real element
%
% INPUT:
%   c:             current cell
%   dim:           dimension of problem (must be 2)
%   numLocalEdges: number of edges in each cell (equal for all cells)
%   nodes:         numNodes by 2 array. The ith row contains the (x,y) 
%                  coordinates of node_i
%   edges:         numEdges by 2 array. The ith row contains pointers into 
%                  nodes, identifying the endpoints of edge_i
%   cells:         numCells by numLocalEdges array.  The ith row contains
%                  pointers into edges, identifying the edges of cell_i, in
%                  counterclockwise order, and their orientations. 
%                  Specifically, suppose cells(i,j)=k.  If k>0, then the 
%                  jth edge of cell_i is edge_k, traced from its first 
%                  endpoint to its second.  If k<0, then the jth edge of 
%                  cell_i is edge_{-k}, traced from its second endpoint to 
%                  its first
% 
% OUTPUT:
%   realNodes: dim by numLocalEdges array containing the nodes of the real 
%              element
%
function realNodes = getRealNodes(c, dim, numLocalEdges, nodes, edges, cells)

realNodes = zeros(dim,numLocalEdges);
for lE=1:numLocalEdges
    edgePtr = abs(cells(lE,c));
    orient = sign(cells(lE,c));
    if orient == 1
        nPtr = edges(1,edgePtr);
    else
        nPtr = edges(2,edgePtr);
    end
    
    realNodes(:,lE) = nodes(:,nPtr);
end
end