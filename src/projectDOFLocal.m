function [uVal, gradUVal, v_f, gradV_f, p_f, gradP_f, v_m, gradV_m, p_m,gradP_m] = projectDOFLocal(u, c, realBasis, realBasisDer, appCtx)

dim = appCtx.dim;
numFields = appCtx.numFields;
numQuadPoints = appCtx.quad.numQuadPoints;

uVal = cell(numQuadPoints,1);
gradUVal = cell(numQuadPoints,1);

v_f = zeros(2,numQuadPoints);
gradV_f = zeros(2,2,numQuadPoints);
p_f = zeros(1,numQuadPoints);
gradP_f = zeros(2,numQuadPoints);
v_m = zeros(2,numQuadPoints);
gradV_m = zeros(2,2,numQuadPoints);
p_m = zeros(1,numQuadPoints);
gradP_m = zeros(2,numQuadPoints);
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
    if f == 1
      v_f(:,q) = uValTmp;
      gradV_f(:,:,q) = gradUValTmp;
    elseif f == 2
      p_f(q) = uValTmp;
      gradP_f(:,q) = gradUValTmp;
    elseif f == 3
      v_m(:,q) = uValTmp;
      gradV_m(:,:,q) = gradUValTmp;
    elseif f == 4
      p_m(q) = uValTmp;
      gradP_m(:,q) = gradUValTmp;
    end
  end
end
