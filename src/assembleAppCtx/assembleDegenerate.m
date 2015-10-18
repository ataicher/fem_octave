function appCtx = assembleDegenerate(phi, Theta, RT0VelocityField, RT0PressureField ,underIntegrate, appCtx)

numFields = appCtx.numFields;

appCtx.degenerate.Theta = Theta;

if isa(phi,'function_handle') && nargin(phi) == 2
    appCtx.degenerate.phi = projectPhiFunction(phi,appCtx);
    appCtx.degenerate.phiFun = phi;
else
    error('myApp:argChk', 'porosity phi is not defined properly. Must be a function handle of two variables')
end

if RT0VelocityField <= numFields && RT0PressureField <= numFields
    appCtx.degenerate.RT0VelocityField = RT0VelocityField;
    appCtx.degenerate.RT0PressureField = RT0PressureField;
else
    error('myApp:argChk', 'RT0VelocityField and RT0PressureField do not have proper values')
end

if strcmp(underIntegrate, 'true')
    appCtx.degenerate.UNDERINTEGRATE = 1;
else
    appCtx.degenerate.UNDERINTEGRATE = 0;
end

end

function phi = projectPhiFunction(phiFun,appCtx)

numCells = appCtx.mesh.numCells;
numEdges = appCtx.mesh.numEdges;
numLocalEdges = appCtx.mesh.numLocalEdges;
cells = appCtx.mesh.cells;
numQuadPoints = appCtx.quad.numQuadPoints;
numEdgeQuadPoints = appCtx.quad.numEdgeQuadPoints;
quadWeights = appCtx.quad.quadWeights;
edgeQuadWeights = appCtx.quad.edgeQuadWeights;
refArea = appCtx.quad.refArea;

% local to global map
phi.cellDOF.c = (1:numCells);
phi.cellDOF.e = abs(cells);

phi.c = zeros(numCells,1);
phi.e = zeros(numEdges,1);
for c=1:numCells
    detJ = appCtx.cellGeometry.detJ{c};
    realArea = detJ*refArea;
    realQuadPoints = projectQuadPoints(c,appCtx);
    
    % cell dof values
    for q=1:numQuadPoints
        phi.c(c) = phi.c(c) + ...
            (1/realArea)*detJ*quadWeights(q)*phiFun(realQuadPoints(1,q),realQuadPoints(2,q));
    end
    
    % edge dof values
    for lEPtr=1:numLocalEdges
        realEdgeQuadPoints = projectEdgeQuadPoints(lEPtr,c,appCtx);
        eQuadWeights = edgeQuadWeights(:,lEPtr);
        ePtr = abs(cells(lEPtr,c));
        phi.e(ePtr) = 0;
        for q=1:numEdgeQuadPoints
            phi.e(ePtr) = phi.e(ePtr) + ...
                (1/sqrt(refArea))*eQuadWeights(q)*phiFun(realEdgeQuadPoints(1,q),realEdgeQuadPoints(2,q));
        end
    end
end
end