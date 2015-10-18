function plotBasisFuncs(appCtx)

numFields = appCtx.numFields;
quadPoints = appCtx.quad.quadPoints;

xGrid = unique(sort(quadPoints(1,:)));
yGrid = unique(sort(quadPoints(2,:)));
[X,Y] = meshgrid(xGrid,yGrid);

for f = 1:numFields
    numComp = appCtx.field(f).numComp;
    
    if numComp == 1
        
        numRows = floor(sqrt(appCtx.field(f).numBasisFuncs));
        if floor(sqrt(appCtx.field(f).numBasisFuncs)) == sqrt(appCtx.field(f).numBasisFuncs)
            numCols = numRows;
        else
            numCols = numRows + ceil((appCtx.field(f).numBasisFuncs-numRows^2)/numRows);
        end
        
    elseif numComp == appCtx.dim
        
        numRows = floor(sqrt(numComp*appCtx.field(f).numBasisFuncs));
        if floor(sqrt(numComp*appCtx.field(f).numBasisFuncs)) == sqrt(numComp*appCtx.field(f).numBasisFuncs)
            numCols = numRows;
        else
            numCols = numRows + ceil((numComp*appCtx.field(f).numBasisFuncs-numRows^2)/numRows);
        end
    end
    
    figure, clf,
    
    for b=1:appCtx.field(f).numBasisFuncs
        for comp = 1:numComp
            fun = reshape(appCtx.field(f).basis(b,comp,:),appCtx.quad.numQuadPoints,1);
            Fun = griddata(appCtx.quad.quadPoints(1,:),appCtx.quad.quadPoints(2,:),fun,X,Y);
            subplot(numRows,numCols,(b-1)*numComp + comp)
            surf(X,Y,Fun)
            xlabel('x')
            ylabel('y')
            zlabel(['basis: ', num2str(b), ' comp: ' num2str(comp)])
        end
    end
    suptitle(['basis functions on the reference element for field ' num2str(f)])
end