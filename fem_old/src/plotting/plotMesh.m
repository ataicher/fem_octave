%% PLOTMESH
% view the mesh with global numbers for nodes, edges, and cells
% INPUT:
%   mesh: mesh structure contatining all mesh information
% OUTPUT:
%   plot: mesh plot showing cell, edge, and node numbering
%
function plotMesh(appCtx, varargin)

numCells = appCtx.mesh.numCells;
numEdges = appCtx.mesh.numEdges;
numNodes = appCtx.mesh.numNodes;
cells = appCtx.mesh.cells;
edges = appCtx.mesh.edges;
nodes = appCtx.mesh.nodes;
numLocalEdges = appCtx.mesh.numLocalEdges;

hold on

% plot edges
for ePtr=1:numEdges
    n1Ptr = edges(1,ePtr);
    n2Ptr = edges(2,ePtr);
    n1 = nodes(:,n1Ptr);
    n2 = nodes(:,n2Ptr);
    plot([n1(1) n2(1)],[n1(2) n2(2)],'k')
end

% viewing region
nodeMaxX = max(nodes(1,:));
nodeMinX = min(nodes(1,:));
dX = norm(nodeMinX-nodeMaxX);
nodeMaxY = max(nodes(2,:));
nodeMinY = min(nodes(2,:));
dY = norm(nodeMinY-nodeMaxY);
spaceX = dX/15;
spaceY = dY/15;
axis([nodeMinX-spaceX nodeMaxX+spaceX nodeMinY-spaceY nodeMaxY+spaceY])

% highlight edge without boundary condition if it exists
if nargin > 1
    ePtr = varargin{1};
    n1Ptr = edges(1,ePtr);
    n2Ptr = edges(2,ePtr);
    n1 = nodes(:,n1Ptr);
    n2 = nodes(:,n2Ptr);
    coords =  (n1 + n2)/2;
    text(coords(1),coords(2),int2str(ePtr),'FontSize',100,'Color','cyan');
    return
end

% cell numbering
for cPtr = 1:numCells
    coords = [0 ; 0];
    for lE=1:numLocalEdges
        ePtr = abs(cells(lE,cPtr));
        n1Ptr = edges(1,ePtr);
        n2Ptr = edges(2,ePtr);
        n1 = nodes(:,n1Ptr);
        n2 = nodes(:,n2Ptr);
        coords = coords + n1 + n2;
    end
    coords = coords/(2*numLocalEdges);
    text(coords(1),coords(2),coords(2),int2str(cPtr),'FontSize',100,'Color','magenta');
end

% edge numbering
for ePtr=1:numEdges
    n1Ptr = edges(1,ePtr);
    n2Ptr = edges(2,ePtr);
    n1 = nodes(:,n1Ptr);
    n2 = nodes(:,n2Ptr);
    coords =  (n1 + n2)/2;
    text(coords(1),coords(2),int2str(ePtr),'FontSize',100,'Color','red');
end

% node numbering
for nPtr=1:numNodes
    coords = [nodes(1,nPtr) ; nodes(2,nPtr)];
    text(coords(1),coords(2),int2str(nPtr),'FontSize',100,'Color','blue');
end

title(['Mesh: ', '{\color{blue}', num2str(numNodes),' nodes}, ', '{\color{red}', num2str(numEdges),...
    ' edges}, ', '{\color{magenta}', num2str(numCells), ' cells}'],'FontSize',80)