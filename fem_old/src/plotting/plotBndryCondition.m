%% PLOTMESH
% view the mesh with global numbers for nodes, edges, and cells
% INPUT:
%   mesh: mesh structure contatining all mesh information
% OUTPUT:
%   plot: mesh plot showing cell, edge, and node numbering
%
function plotBndryCondition(appCtx, varargin)

numFields = appCtx.numFields;
numBndryEdges = appCtx.mesh.numBndryEdges;
bndryEdges = appCtx.mesh.bndryEdges;
edges = appCtx.mesh.edges;
nodes = appCtx.mesh.nodes;

fieldBC = [];
for f=1:numFields
    if appCtx.field(f).BNDRYCONDTION
        fieldBC = [fieldBC f];
    end
end
numFieldBC = length(fieldBC);

iter = 1;
for f = fieldBC
    name = appCtx.field(f).name;
    subplot(numFieldBC,1,iter)
    hold on
    iter = iter + 1;
    numBC = appCtx.field(f).numBC;
    edgeBC = appCtx.field(f).edgeBC;
    cc = hsv(numBC);
    h = zeros(numBC,1);
    % plot edges
    for bEPtr=1:numBndryEdges
        ePtr = bndryEdges(bEPtr);
        n1Ptr = edges(1,ePtr);
        n2Ptr = edges(2,ePtr);
        n1 = nodes(:,n1Ptr);
        n2 = nodes(:,n2Ptr);
        i = edgeBC(bEPtr);
        h(i) = plot([n1(1) n2(1)],[n1(2) n2(2)],'color',cc(i,:),'LineWidth', 3);
    end
    
    for i=1:numBC
        legendStr{i} = ['bndry' num2str(i)];
    end
    legend(h, legendStr)
    clearvars legendStr;
    
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
    
    title(name,'FontSize', 15)
    
end

suptitle('Boundary Coniditons')
legH = legend;
if ~isempty(legH)
    axes(legH);
end
