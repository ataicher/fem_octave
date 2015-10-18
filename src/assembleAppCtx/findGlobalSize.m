% Get total number of degrees of freedom
% INPUT:
% OUTPUT:
%   size - global number of degrees of freedom
function [globalSize, fieldDOF] = findGlobalSize(appCtx)

numFields = appCtx.numFields;
numNodes = appCtx.mesh.numNodes;
numEdges = appCtx.mesh.numEdges;
numCells = appCtx.mesh.numCells;

for f=1:numFields
    numNodeDOF = appCtx.field(f).numNodeDOF;
    numEdgeDOF = appCtx.field(f).numEdgeDOF;
    numFaceDOF = appCtx.field(f).numFaceDOF;
    
    fieldDOF(f) = numNodes*numNodeDOF + numEdges*numEdgeDOF + numCells*numFaceDOF;
end

globalSize = sum(fieldDOF);
