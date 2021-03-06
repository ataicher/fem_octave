% UNIFORMQUADMESH
% Create a uniform quadrelateral mesh.
% INPUT:
%   Lx - width of domain
%   Ly - height of domain
%   nx - number of cells in the x-direction
%   ny - number of cells in the y-direction
% OUTPUT:
%   cells: numCells by numLocalEdges array.  The ith row contains
%          pointers into edges, identifying the edges of cell_i, in 
%          counterclockwise order, and their orientations. Specifically, 
%          suppose cells(i,j)=k.  If k>0, then the jth edge of cell_i is 
%          edge_k, traced from its first endpoint to its second.  If k<0, 
%          then the jth edge of cell_i is edge_{-k}, traced from its second
%          endpoint to its first.
%
%   edges: numEdges by 2 array. The ith row contains pointers into nodes, 
%          identifying the endpoints of edge_i.
%
%   nodes: numNodes by 2 array.  The ith row contains the (x,y) coordinates
%          of node_i.
%
function [nodes, edges, cells] = quadMesh(x,y,varargin)

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
% horizontal edges
for j=1:ny+1
    for i=1:nx
        ePtr = (j-1)*nx + i;
        n1Ptr = (j-1)*(nx+1) + i;
        n2Ptr = (j-1)*(nx+1) + i + 1;
        edges(:,ePtr) = [n1Ptr ; n2Ptr];
% check if edge is connected to a removed node
        if sum((nRemovePtrs == n1Ptr) + (nRemovePtrs == n2Ptr)) > 0
            eRemovePtrs = [eRemovePtrs ePtr];
        end
    end
end

% vertical edges
for i=1:nx+1
    for j=1:ny
        ePtr = (ny+1)*nx + (i-1)*ny + j;
        n1Ptr = (j-1)*(nx+1) + i;
        n2Ptr = j*(nx+1) + i;
        edges(:,ePtr) = [n1Ptr ; n2Ptr];
% check if edge is connected to a removed node
        if sum((nRemovePtrs == n1Ptr) + (nRemovePtrs == n2Ptr)) > 0
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
        if sum((eRemovePtrs == edgeBPtr) + (eRemovePtrs == edgeRPtr) + (eRemovePtrs == edgeTPtr) + (eRemovePtrs == edgeLPtr)) > 0
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
