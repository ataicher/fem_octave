function plotFunDer(xPoints, yPoints, funVal, funDerVal, appCtx)
dim = appCtx.dim;
variables = appCtx.variables;

numFields = variables.numFields;
numCompTotal = variables.numCompTotal;

xGrid = unique(sort(xPoints));
yGrid = unique(sort(yPoints));
[X,Y] = meshgrid(xGrid,yGrid);
numPoints = length(xPoints);

compInd = 0;
for f=1:numFields
    numComp = variables.field{f}.numComp;
    name = variables.field{f}.name;
    for comp=1:numComp
        FuncVal = griddata(xPoints,yPoints,funVal.field{f}(comp,:)',X,Y);
        compInd = compInd + 1;
        subplot(numCompTotal, dim+1 ,(compInd-1)*3 + 1)
        surf(X,Y,FuncVal),
        xlabel('x')
        ylabel('y')
        if numComp > 1
            title([name, ' comp ',num2str(comp)],'FontSize', 15)
        else
            title(name,'FontSize', 15)
        end
        for d=1:dim
            fDer = reshape(funDerVal.(fName)(d,comp,:),numPoints,1);
            FuncDerVal = griddata(xPoints,yPoints,reshape(funDerVal.(fName)(d,comp,:),numPoints,1),X,Y);
            subplot(numCompTotal, dim+1 ,(compInd-1)*3+d+1)
            surf(X,Y,FuncDerVal),
            xlabel('x')
            ylabel('y')
            if numComp > 1
                title([name, ' comp ',num2str(comp), ' deriv ' , num2str(d)],'FontSize', 15)
            else
                title([name, ' deriv ' , num2str(d)],'FontSize', 15)
            end
        end
    end
end