%% ASSEMBLE LOCAL TO GLOBAL MAP
% Create local to global mapping matrix
%
% INPUT:
%  appCtx: struct containing problem information
%
% OUTPUT:
%   LG: local to global mapping. LG.field{f}(localDof,cell) = global degree
%       of freedom of the local DOF for the cell
%
function appCtx = assembleLocalToGlobalMap(appCtx)

numFields = appCtx.numFields;
edges = appCtx.mesh.edges;
cells = appCtx.mesh.cells;
numNodes = appCtx.mesh.numNodes;
numEdges = appCtx.mesh.numEdges;
numCells = appCtx.mesh.numCells;
numLocalEdges = appCtx.mesh.numLocalEdges;

fInd = 0;
for f=1:numFields
    numNodeDOF = appCtx.field(f).numNodeDOF;
    numEdgeDOF = appCtx.field(f).numEdgeDOF;
    numFaceDOF = appCtx.field(f).numFaceDOF;
    numBasisFuncs = appCtx.field(f).numBasisFuncs;
    
    appCtx.field(f).startDOF = fInd + 1;
    cellDOF = zeros(numBasisFuncs,numCells);
    for c=1:numCells
        cellDOFLocal = zeros(numBasisFuncs,1);
        ind = 1;
        
        if numNodeDOF > 0
            for lEPtr=1:numLocalEdges
                ePtr = abs(cells(lEPtr,c));
                orient = sign(cells(lEPtr,c));
                if orient == 1
                    nPtr = edges(1,ePtr);
                else
                    nPtr = edges(2,ePtr);
                end
                for nDOF=1:numNodeDOF
                    cellDOFLocal(ind) = fInd + (nPtr-1)*numNodeDOF + nDOF;
                    ind = ind + 1;
                end
            end
        end
        
        if numEdgeDOF > 0
            for lEPtr=1:numLocalEdges
                ePtr = abs(cells(lEPtr,c));
                for eDOF=1:numEdgeDOF
                    cellDOFLocal(ind) = fInd + numNodes*numNodeDOF + (ePtr-1)*numEdgeDOF + eDOF;
                    ind = ind + 1;
                end
            end
        end
        
        if numFaceDOF > 0
            for fDOF=1:numFaceDOF
                cellDOFLocal(ind) = fInd + numNodes*numNodeDOF + numEdges*numEdgeDOF + (c-1)*numFaceDOF + fDOF;
                ind = ind + 1;
            end
        end
        
        cellDOF(:,c) = cellDOFLocal; 
    end
    
    fInd = fInd + numNodes*numNodeDOF + numEdges*numEdgeDOF + numCells*numFaceDOF;
    appCtx.field(f).endDOF = fInd;
    appCtx.field(f).cellDOF = cellDOF;
end
