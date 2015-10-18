function condJacVals = compactingColumnZeroPorosityLimit(SETPHI, nx, phiVec, FORMULATION)

% get problem info
appCtx.storage.phiVec = phiVec;
appCtx.storage.SETPHI = SETPHI;
appCtx.storage.FORMULATION = FORMULATION;
appCtx.nx = nx;
compactingColumn;
assembleAppCtx;
EXISTEXACTSOL = appCtx.EXISTEXACTSOL;
phiVec = appCtx.storage.phiVec;
SETPHI = appCtx.storage.SETPHI;
FORMULATION = appCtx.storage.FORMULATION;

% get fluid pressure and relative velocity fields
fFluidPres = 2;
fRelVel = 1;

figure, clf,

lenPhiVec = length(phiVec);

% set up legend and plot coloring
cc = hsv(lenPhiVec+(EXISTEXACTSOL>0));
h = zeros(lenPhiVec+(EXISTEXACTSOL>0),1);
legendStr = cell(lenPhiVec+(EXISTEXACTSOL>0),1);
for i = 1:lengthPhiVec+(EXISTEXACTSOL>0)
    if i==1
        MarkerStyle = 'o-';
    elseif i==2
        MarkerStyle = '*-';
    elseif i==3
        MarkerStyle = 'x-';
    else
        MarkerStyle = 's-';
    end
end

condJacVals = [];

for i = 1:lenPhiVec
    
    % run driver for current value of porosity
    [u, ~, ~, appCtx] = exec(@compactingColumn,0,0,40,1,SETPHI,phiVec(i),phiVec(1),FORMULATION);
    
    condJacVals = [condJacVals appCtx.condJac];
    
    % get values of computed solution at quadrature points
    [xPoints, yPoints, uVal] = projectDOF(u, appCtx);
    xGrid = unique(sort(xPoints));
    yGrid = unique(sort(yPoints));
    [X,Y] = meshgrid(xGrid,yGrid);
    
    % plot fluid pressure or scaled fluid pressure
    funPlot = griddata(xPoints,yPoints,uVal.field{fFluidPres}(1,:)',X,Y);
    funPlot = funPlot(1,:);
    
    subplot(1,5,[4 5]), hold on
    h(i) = line_fewer_markers(xGrid, funPlot, 10, MarkerStyle(i),'MarkerSize', 8,'color',cc(i,:),'LineWidth', 3);
    view(90,90)
    xlabel('depth','FontSize', 15)
    title(appCtx.field(fFluidPres).name,'FontSize', 17)
    
    % plot fluid pressure if degenerate code
    if appCtx.DEGENERATE
        pDOF = u;
        for c=1:appCtx.mesh.numCells
            cellDOF = appCtx.field(fFluidPres).cellDOF(c);
            uLocal = u(cellDOF);
            if appCtx.degenerate.phi.c(c) > 1e-10
                pDOF(cellDOF) = (1/sqrt(appCtx.degenerate.phi.c(c)))*uLocal;
            else
                pDOF(cellDOF) = 0;
            end
        end
        [xPoints, yPoints, pVal] = projectDOF(pDOF, appCtx);
        funPlot = griddata(xPoints,yPoints,pVal.field{fFluidPres}(1,:)',X,Y);
        funPlot = funPlot(1,:);
        subplot(1,5,[4 5]), hold on
        plot(xGrid,funPlot,'color',cc(i,:),'LineWidth', 3);
        view(90,90)
        title('fluid pressure','FontSize', 17)
    end
    
    % plot relative velocity
    funPlot = griddata(xPoints,yPoints,uVal.field{fRelVel}(1,:)',X,Y);
    funPlot = funPlot(1,:);
    subplot(1,5,[2 3]), hold on
    line_fewer_markers(xGrid, funPlot, 10, MarkerStyle,'MarkerSize', 8,'color',cc(i,:),'LineWidth', 3);
    %plot(xGrid,funPlot,'color',cc(i,:),'LineWidth', 3);
    view(90,90)
    title('relative velocity','FontSize', 17)
    
    % plot porosity
    phiVal = zeros(1,length(xGrid));
    for j = 1:length(xGrid)
        phiVal(j) = appCtx.phi(xGrid(j),yGrid(1));
    end
    subplot(1,5,1), hold on
    plot(xGrid, phiVal,'color',cc(i,:), 'LineWidth', 3)
    view(90,90)
    title('porosity','FontSize', 17)
    
    % update legend
    legendStr{i} = num2str(phiVec(i));
end

if EXISTEXACTSOL
    
    % get values of exact solution at quadrature points
    [xPoints, yPoints, uExactVal] = evalExactSolution(appCtx);
    
    % plot exact fluid pressure
    subplot(1,5,[4 5]), hold on
    funPlot = griddata(xPoints,yPoints,uExactVal.field{fFluidPres}(1,:)',X,Y);
    funPlot = funPlot(1,:);
    h(lenPhiVec+1) = plot(xGrid, funPlot,'color',cc(lenPhiVec+1,:),'LineWidth', 3,'LineStyle','--');
    legendStr{lenPhiVec+1} = 'exact solution';
    
    % plot exact relative velocity
    subplot(1,5,[2 3]), hold on
    funPlot = griddata(xPoints,yPoints,uExactVal.field{fRelVel}(1,:)',X,Y);
    funPlot = funPlot(1,:);
    plot(xGrid, funPlot,'color',cc(lenPhiVec+1,:),'LineWidth', 3,'LineStyle','--');
end

legend(h, legendStr);