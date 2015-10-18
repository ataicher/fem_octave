function realEdgeQuadPoints = projectEdgeQuadPoints(lEPtr,c,appCtx)

dim = appCtx.dim;
edgeQuadPoints = appCtx.quad.edgeQuadPoints(:,:,lEPtr);
numEdgeQuadPoints = appCtx.quad.numEdgeQuadPoints;
v_0Ref = appCtx.quad.v_0Ref;
v_0 = appCtx.cellGeometry.v_0{c};
J = appCtx.cellGeometry.J{c};

realEdgeQuadPoints = zeros(dim,numEdgeQuadPoints);
for q = 1:numEdgeQuadPoints
    realEdgeQuadPoints(:,q) = J*(edgeQuadPoints(:,q) - v_0Ref) + v_0;
end