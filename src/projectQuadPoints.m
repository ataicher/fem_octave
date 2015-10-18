function realQuadPoints = projectQuadPoints(c,appCtx)

dim = appCtx.dim;
quadPoints = appCtx.quad.quadPoints;
numQuadPoints = appCtx.quad.numQuadPoints;
v_0Ref = appCtx.quad.v_0Ref;
v_0 = appCtx.cellGeometry.v_0{c};
J = appCtx.cellGeometry.J{c};

realQuadPoints = zeros(dim,numQuadPoints);
for q = 1:numQuadPoints
    realQuadPoints(:,q) = J*(quadPoints(:,q) - v_0Ref) + v_0;
end
