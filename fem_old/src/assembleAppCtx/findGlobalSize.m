% Get total number of degrees of freedom
% INPUT:
% OUTPUT:
%   size - global number of degrees of freedom
function globalSize = findGlobalSize(appCtx)

numFields = appCtx.numFields;
numNodes = appCtx.mesh.numNodes;
numEdges = appCtx.mesh.numEdges;
numCells = appCtx.mesh.numCells;

globalSize = 0;
for f=1:numFields
    numNodeDOF = appCtx.field(f).numNodeDOF;
    numEdgeDOF = appCtx.field(f).numEdgeDOF;
    numFaceDOF = appCtx.field(f).numFaceDOF;
    
    globalSize = globalSize + numNodes*numNodeDOF + numEdges*numEdgeDOF + numCells*numFaceDOF;
end
