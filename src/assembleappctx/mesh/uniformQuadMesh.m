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
function [nodes edges cells] = uniformQuadMesh(nx,ny,Lx,Ly)

if nargin==3
   Ly=Lx;
elseif nargin<3
   Lx=[-1 1];
   Ly=[-1 1];
end
if nargin<2
   ny=nx;
end

% NODES
numNodes = (nx+1)*(ny+1);
nodes = zeros(2,numNodes);
for j=1:ny+1
    for i=1:nx+1
        nPtr = (j-1)*(nx+1)+i;
        xCoord = (i-1)*(Lx(2)-Lx(1))/nx + Lx(1);
        yCoord = (j-1)*(Ly(2)-Ly(1))/ny + Ly(1);
        nodes(:,nPtr) = [xCoord ; yCoord];
    end
end

% EDGES
numEdges = (ny+1)*nx+ny*(nx+1);
edges = zeros(2,numEdges);

% horizontal edges
for j=1:ny+1
    for i=1:nx
        ePtr = (j-1)*nx + i;
        n1Ptr = (j-1)*(nx+1) + i;
        n2Ptr = (j-1)*(nx+1) + i + 1;
        edges(:,ePtr) = [n1Ptr ; n2Ptr];
    end
end

% vertical edges
for i=1:nx+1
    for j=1:ny
        ePtr = (ny+1)*nx + (i-1)*ny + j;
        n1Ind = (j-1)*(nx+1) + i;
        n2Ind = j*(nx+1) + i;
        edges(:,ePtr) = [n1Ind ; n2Ind];
    end
end

% CELLS
numCells = nx*ny;
cells = zeros(4,numCells);
for j=1:ny
    for i=1:nx
        edgeBPtr = (j-1)*nx + i;
        edgeRPtr = nx*(ny+1) + i*ny + j;
        edgeTPtr = j*nx + i;
        edgeLPtr = nx*(ny+1) + (i-1)*ny + j;
        cells(:,(j-1)*nx+i) = [edgeBPtr edgeRPtr -edgeTPtr -edgeLPtr]';
    end
end