function plotCompactingColumn(u, appCtx)

numFields = appCtx.numFields;
phi = appCtx.phi;
EXISTEXACTSOL = appCtx.EXISTEXACTSOL;
DEGENERATE = appCtx.DEGENERATE;

    % get values of computed and exact (if exists) solution at quadrature points
[xPoints, yPoints, uVal] = projectDOF(u, appCtx);
if EXISTEXACTSOL
    [xPoints, yPoints, uExactVal] = evalExactSolution(appCtx);
end

figure, clf,

xGrid = unique(sort(xPoints));
yGrid = unique(sort(yPoints));
[X,Y] = meshgrid(xGrid,yGrid);

% plot values for each field
for f=1:numFields 
    name = appCtx.field(f).name;

        funPlot = griddata(xPoints,yPoints,uVal.field{f}(1,:)',X,Y);
        funPlot = funPlot(1,:);
        
        subplot(1,numFields+3,f)
        plot(xGrid, funPlot, 'LineWidth', 4,'color','black')
        
        if EXISTEXACTSOL
            hold on
            funPlot = griddata(xPoints,yPoints,uExactVal.field{f}(1,:)',X,Y);
            funPlot = funPlot(1,:);
            plot(xGrid, funPlot, 'LineWidth', 4,'color','red','LineStyle','--')
            
        end
                     ax = gca;
            set(ax,'Box','on','FontSize',15);       
        view(90,90)
        if f == 1
            xlabel('depth','FontSize', 20)
            if EXISTEXACTSOL
                legend('computed','exact','Location','NorthWest')
            else
                legend('computed','Location','NorthWest')
            end
        else
            set(ax,'XTick',[])
        end
        title(name,'FontSize', 21)
        
end

% plot scaled fluid or fluid pressure
if 1
if DEGENERATE
    pDOF = u;
    for c=1:appCtx.mesh.numCells
        cellDOF = appCtx.field(2).cellDOF(c);
        uLocal = u(cellDOF);
        if appCtx.degenerate.phi.c(c) > 1e-10
            pDOF(cellDOF) = (1/sqrt(appCtx.degenerate.phi.c(c)))*uLocal;
        end
    end
    [xPoints, yPoints, pVal] = projectDOF(pDOF, appCtx);
    
    funPlot = griddata(xPoints,yPoints,pVal.field{2}(1,:)',X,Y);
    funPlot = funPlot(1,:);
    %     for i = 1:length(xGrid)
    %         if phi(xGrid(i),yGrid(1)) > 1e-10
    %             qVal(i) = (1/sqrt(phi(xGrid(i),yGrid(1))))*funPlot(i);
    %         else
    %             qVal(i) = 0;
    %         end
    %     end
else
    funPlot = griddata(xPoints,yPoints,uVal.field{2}(1,:)',X,Y);
    funPlot = funPlot(1,:);
    for i = 1:length(xGrid)
        funPlot(i) = sqrt(phi(xGrid(i),yGrid(1)))*funPlot(i);
    end
end
subplot(1,numFields+3,numFields+2)
plot(xGrid, funPlot, 'LineWidth', 4,'color','black')
view(90,90)
ax = gca;
set(ax,'Box','on','FontSize',15,'XTick',[]);
if DEGENERATE
    title('fluid pressure','FontSize', 21)
else
    title('scaled pressure','FontSize', 21)
end
end

% plot porosity
phiVal = zeros(1,length(xGrid));
for i = 1:length(xGrid)
    phiVal(i) = phi(xGrid(i),yGrid(1));
end
subplot(1,numFields+3,numFields+1)
plot(xGrid, phiVal, 'LineWidth', 4)
view(90,90)
ax = gca;
set(ax,'Box','on','FontSize',15,'XTick',[]);
title('porosity','FontSize', 21)