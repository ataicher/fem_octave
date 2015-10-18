function [uVal, gradUVal] = projectDOFLocal(u, c, realBasis, realBasisDer, appCtx)

dim = appCtx.dim;
numFields = appCtx.numFields;
numQuadPoints = appCtx.quad.numQuadPoints;

uVal = cell(numQuadPoints,1);
gradUVal = cell(numQuadPoints,1);
for f=1:numFields
    numComp = appCtx.field(f).numComp;
    rBasis = realBasis.field{f};
    rBasisDer = realBasisDer.field{f};
    cellDOF = appCtx.field(f).cellDOF(:,c);
    
    uLocal = u(cellDOF);
    
    for q=1:numQuadPoints
        uValTmp = zeros(numComp,1);
        gradUValTmp = zeros(numComp,dim);
        for comp=1:numComp
                uValTmp(comp) = uValTmp(comp) + uLocal'*rBasis(:,comp,q);
                for d=1:dim
                    gradUValTmp(comp,d) = gradUValTmp(comp,d) + uLocal'*rBasisDer(:,d,comp,q);
                end
        end
        
        uVal{q}.field{f} = uValTmp;
        gradUVal{q}.field{f} = gradUValTmp;
    end
end