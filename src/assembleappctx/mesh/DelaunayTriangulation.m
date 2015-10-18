function [nodes, edges, cells] = DelaunayTriangulation(X,Y)

dt = DelaunayTri(X,Y);
nodes = dt.X';
edges = dt.edges';
cellNodes = dt.Triangulation';
numNodes = size(nodes,2);
numEdges = size(edges,2);
numCells = size(cellNodes,2);

cells = zeros(3,numCells);
for c=1:numCells
    for e=1:numEdges
        if (edges(1,e) == cellNodes(1,c)) && (edges(2,e) == cellNodes(2,c))
            cells(1,c) = e;
        end
        if (edges(1,e) == cellNodes(2,c)) && (edges(2,e) == cellNodes(1,c))
            cells(1,c) = -e;
        end
        
        if (edges(1,e) == cellNodes(2,c)) && (edges(2,e) == cellNodes(3,c))
            cells(2,c) = e;
        end
        if (edges(1,e) == cellNodes(3,c)) && (edges(2,e) == cellNodes(2,c))
            cells(2,c) = -e;
        end
          
        if (edges(1,e) == cellNodes(3,c)) && (edges(2,e) == cellNodes(1,c))
            cells(3,c) = e;
        end
        if (edges(1,e) == cellNodes(1,c)) && (edges(2,e) == cellNodes(3,c))
            cells(3,c) = -e;
        end
    end
end