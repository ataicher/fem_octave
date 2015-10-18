%% ASSEMBLE MESH-----------------------------------------------------------
% Assemble mesh information
% INPUT:
%   nodes: numNodes by 2 array. The ith row contains the (x,y) coordinates
%          of node_i
%   edges: numEdges by 2 array. The ith row contains pointers into nodes,
%          identifying the endpoints of edge_i
%   cells: numCells by numLocalEdges array.  The ith row contains
%          pointers into edges, identifying the edges of cell_i, in
%          counterclockwise order, and their orientations. Specifically,
%          suppose cells(i,j)=k.  If k>0, then the jth edge of cell_i is
%          edge_k, traced from its first endpoint to its second.  If k<0,
%          then the jth edge of cell_i is edge_{-k}, traced from its second
%          endpoint to its first
%   viewMesh: 'true' if you want to view the mesh
%
% OUTPUT: mesh struct containing
%   type:          type of mesh (triangle or quad)
%   nodes:         same as input
%   numNodes:      the number of nodes
%   edges:         same as input
%   numEdges:      the number of edges
%   cells:         same as input
%   numCells:      the number of cells
%   numLocalEdges: number of edges in each cell (equal for all cells)
%   edgeCells:     2 by numEdges array. The ith row contains pointers into 
%                  cells, identifying the cells on either side ofedge_i. If
%                  edge_i is a boundary edge, the secondpointer is 0
%   edgeNormals:    dim by numEdges array. The ith column contains the normal
%                  vector for the ith edge
%   bndyEdges:     numBndryEdges long array. The ith entry is the index of 
%                  in edges of the ith boundary edge
%   numBndryEdges: the number of boundary edges
%   viewMesh:  'true' if you want to view the mesh
%
function mesh = assembleMesh(x,y,loc)

dim = 2;

% build mesh
[nodes, edges, cells] = quadMesh(x,y,loc);
numNodes = size(nodes,2);
numEdges = size(edges,2);
numCells = size(cells,2);
numLocalEdges = size(cells,1);

%% INPUT TESTS-------------------------------------------------------------

% nodes size
if size(nodes,1) ~= dim
        error('myApp:argChk', 'nodes must be dimx(numNodes)')
end

% edges size
if size(edges,1) ~= dim
        error('myApp:argChk', 'edges must be dimx(numEdges)')
end

% cells size
if numLocalEdges ~= 4
        error('myApp:argChk', 'cells must be 4x(numEdges), only supports rectangular elements.')
end

% compatitibily between edges and cells
if norm(sort(cells(:)) - unique(cells(:))) > 1e-10 || length(unique(abs(cells(:)))) ~= numEdges
    error('myApp:argChk', 'incompatible cells and edges arrays')
end

% no edge can be represented twice with the same orientation
if length(unique(cells(:))) ~= length(cells(:))
    error('myApp:argChk', 'incompatible cells array: no edge can be represented twice with the same orientation')
end

% all edges must be represented
if sort(unique(abs(cells(:)))) ~= (1:numEdges)'
    error('myApp:argChk', 'incompatible cells and edges arrays: all edges must be represented')
end

% each cell must have unique edges
for c =1:size(cells,2)
    if length(unique(abs(cells(:,c)))) < numLocalEdges
        error('myApp:argChk', 'incompatible cells array: each cell must have unique edges')
    end
end

% check edges have positive length
for ePtr = 1:numEdges
    n1Ptr = edges(1,ePtr);
    n2Ptr = edges(2,ePtr);
    if nodes(:,n1Ptr) == nodes(:,n2Ptr)
        error('myApp:argChk', 'edge %d has no length', ePtr)
    end
end

% check cell edges are connected
for c = 1:numCells
    for lEPtr=1:numLocalEdges;
        e1Ptr = abs(cells(lEPtr,c));
        o1 = sign(cells(lEPtr,c));
        e2Ptr = abs(cells(mod(lEPtr,4)+1,c));
        o2 = sign(cells(mod(lEPtr,4)+1,c));
        if o1 == 1
            n1Ptr = edges(2,e1Ptr);
        else
            n1Ptr = edges(1,e1Ptr);
        end
        if o2 == 1
            n2Ptr = edges(1,e2Ptr);
        else
            n2Ptr = edges(2,e2Ptr);
        end
        
        if nodes(:,n1Ptr) ~= nodes(:,n2Ptr)
            error('myApp:argChk', 'incompatible cells, edges, and nodes arrays: cell edges are not connected')
        end  
    end
end
    
%% ASSEMBLE MESH STRUCT----------------------------------------------------
mesh.nodes = nodes;
mesh.numNodes = numNodes;
mesh.edges = edges;
mesh.numEdges = numEdges;
mesh.cells = cells;
mesh.numCells = numCells;
mesh.numLocalEdges = numLocalEdges;
mesh.edgeCells = edgeToCellMap(cells, numCells, numEdges, numLocalEdges);
mesh.edgeNormals = constructEdgeNormals(nodes, edges, numEdges);
mesh.bndryEdges = find(mesh.edgeCells(2,:)==0);
mesh.numBndryEdges = length(mesh.bndryEdges);
end

%% QUADMESH----------------------------------------------------------------
% assemble rectangular mesh
function [nodes, edges, cells] = quadMesh(x,y,loc)

nx = length(x)-1;
ny = length(y)-1;

nRemovePtrs = []; eRemovePtrs = []; cReomvePtrs = [];

% NODES
numNodes = (nx+1)*(ny+1);
nodes = zeros(2,numNodes);
for j=1:ny+1
    for i=1:nx+1
        nPtr = (j-1)*(nx+1) + i;
        nodes(:,nPtr) = [x(i) ; y(j)];
	% check if node is inside domain
        if (loc(x(i),y(j)) <= 0)
            nRemovePtrs = [nRemovePtrs nPtr];
        end
    end
end

% EDGES
numEdges = (ny+1)*nx + (nx+1)*ny;
edges = zeros(2,numEdges);

% vertical edges

    for j=1:ny
      for i=1:nx+1
        ePtr = (j-1)*(nx+1) + i;
        n1Ptr = (j-1)*(nx+1) + i;
        n2Ptr = j*(nx+1) + i;
        edges(:,ePtr) = [n1Ptr ; n2Ptr];
	% check if edge is connected to a removed node
        if any(nRemovePtrs == n1Ptr) | any(nRemovePtrs == n2Ptr)
            eRemovePtrs = [eRemovePtrs ePtr];
        end
    end
end

% horizontal edges
for j=1:ny+1
    for i=1:nx
        ePtr = (nx+1)*ny + (j-1)*nx + i;
        n1Ptr = (j-1)*(nx+1) + i;
        n2Ptr = (j-1)*(nx+1) + i + 1;
        edges(:,ePtr) = [n1Ptr ; n2Ptr];
% check if edge is connected to a removed node
        if any(nRemovePtrs == n1Ptr) | any(nRemovePtrs == n2Ptr)
            eRemovePtrs = [eRemovePtrs ePtr];
        end
    end
end

% CELLS
numCells = nx*ny;
cells = zeros(4,numCells);
for j=1:ny
    for i=1:nx
        cPtr = (j-1)*nx + i;
        edgeBPtr = (j-1)*nx + i;
        edgeRPtr = (ny+1)*nx + i*ny + j;
        edgeTPtr = j*nx + i;
        edgeLPtr = (ny+1)*nx + (i-1)*ny + j;
        cells(:,cPtr) = [edgeBPtr edgeRPtr -edgeTPtr -edgeLPtr]';
	% check if cell is connected to a removed edge
        if any(eRemovePtrs == edgeBPtr) | any(eRemovePtrs == edgeRPtr) | any(eRemovePtrs == edgeTPtr) | any(eRemovePtrs == edgeLPtr)
            cReomvePtrs = [cReomvePtrs cPtr];
        end
    end
end

% remove nodes
nodes(:,nRemovePtrs) = []; numNodes = size(nodes,2);

% remove edges
edges(:,eRemovePtrs) = []; numEdges = size(edges,2);
for ePtr = 1:numEdges
    edges(1,ePtr) = edges(1,ePtr) - sum(edges(1,ePtr)>nRemovePtrs);
    edges(2,ePtr) = edges(2,ePtr) - sum(edges(2,ePtr)>nRemovePtrs);
end

% remove cells
cells(:,cReomvePtrs) = []; numCells = size(cells,2);
for cPtr = 1:numCells
    cells(1,cPtr) = cells(1,cPtr) - sign(cells(1,cPtr))*sum(abs(cells(1,cPtr))>eRemovePtrs);
    cells(2,cPtr) = cells(2,cPtr) - sign(cells(2,cPtr))*sum(abs(cells(2,cPtr))>eRemovePtrs);
    cells(3,cPtr) = cells(3,cPtr) - sign(cells(3,cPtr))*sum(abs(cells(3,cPtr))>eRemovePtrs);
    cells(4,cPtr) = cells(4,cPtr) - sign(cells(4,cPtr))*sum(abs(cells(4,cPtr))>eRemovePtrs);
end
end

%% EDGETOCELLS-------------------------------------------------------------
% Construct edgeCells array
% INPUT: 
%   cells:         numCells by numLocalEdges array.  The ith row contains 
%                  pointers into edges, identifying the edges of cell_i, in
%                  counterclockwise order, and their orientations. 
%                  Specifically, suppose cells(i,j)=k.  If k>0, then the 
%                  jth edge of cell_i is edge_k, traced from its first 
%                  endpoint to its second.  If k<0, then the jth edge of 
%                  cell_i is edge_{-k}, traced from its second endpoint to 
%                  its first
%   numCells:      number of cells
%   numEdges:      number of edges
%   numLocalEdges: number of edges in each cell (equal for all cells)
%
% OUTPUT: 
%   edgeCells: 2 by numEdges array. The ith row contains pointers into 
%              cells, identifying the cells on either side ofedge_i. If
%              edge_i is a boundary edge, the secondpointer is 0
%
function edgeCells = edgeToCellMap(cells, numCells, numEdges, numLocalEdges)

edgeCells = zeros(2,numEdges);
for c=1:numCells
    for lEPtr=1:numLocalEdges
        ePtr = abs(cells(lEPtr,c));
        edgeCells((edgeCells(1,ePtr)>0)+1,ePtr) = c;
    end
end

end

%% CONSTRUCTEDGENORMALS----------------------------------------------------
% Construct edgeNormals array
% INPUT: 
%   nodes:    2 by numNodes array. The ith row contains the (x,y) 
%             coordinates of node_i
%   edges:    2 by numEdges  array. The ith row contains pointers into 
%             nodes, identifying the endpoints of edge_i
%   numEdges: number of edges
%
% OUTPUT: 
%   edgeNormals: dim by numEdges array. The ith column contains the normal
%                 vector for the ith edge
%
function edgeNormals = constructEdgeNormals(nodes, edges, numEdges)

edgeNormals = zeros(2,numEdges);
for ePtr=1:numEdges
    n1Ptr = edges(1,ePtr);
    n2Ptr = edges(2,ePtr);
    e = nodes(:,n2Ptr) - nodes(:,n1Ptr);
    edgeNormals(:,ePtr) = [0 1 ; -1 0]*e/norm(e);
end
end
