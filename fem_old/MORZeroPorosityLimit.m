function MORZeroPorosityLimit(phiVec, nx, ny, SETPHI)

switch SETPHI
    case 1
        fprintf('\nconstant prosity: phi = B\n')
    case 2
        fprintf('\nMOR like porosity')
        fprintf('      { phi_0                 y >= y1  ')
        fprintf('phi = { B + A*(x-C)^2   C < x <= (D-B)/A + C\n')
        fprintf('      { D                   x >  (D-B)/A + C\n')
end
fprintf('\n')

% get problem info
MOR;
appCtx.dim = 2;
appCtx.field = field;
appCtx = assembleVariables(appCtx);
appCtx = assembleExactSolution(appCtx);
EXISTEXACTSOL = appCtx.EXISTEXACTSOL;
numFields = appCtx.numFields;

% get fluid pressure and relative velocity fields
for f = 1:numFields
    name = appCtx.field(f).name;
    if strcmp(name,'fluid pressure potential')
        fFluidPres = f;
    end
    if strcmp(name,'relative velocity')
        fRelVel = f;
    end
end

figure, clf,

lenPhiVec = length(phiVec);

% % set up legend and plot coloring
% if EXISTEXACTSOL
%     cc = hsv(lenPhiVec+1);
%     h = zeros(lenPhiVec+1,1);
%     legendStr = cell(lenPhiVec+1,1);
% else
%     cc = hsv(lenPhiVec);
%     h = zeros(lenPhiVec,1);
%     legendStr = cell(lenPhiVec,1);
% end

for i = 1:lenPhiVec
    
    % run driver for this refinement
    [u, norms, ~, appCtx] = driver(@compactingColumn,0,0,nx,ny,SETPHI,phiVec(i));
    
    % get L2 norm of fluid pressure potential
    normL2 = norms(fFluidPres).L2;
    
    % get scaled L2 norm with porosity
    scaledNorm = L2ScaledNorm(u, fFluidPres, appCtx);
    
    % get problem parameters
    phi = appCtx.phi;
    dxPhi = appCtx.dxPhi;
    n = appCtx.n;
    
    % get Linfty, L1, L2 norms for porosity function
    %   fun = phi^{(n-3)/2}*grad phi
    fun = @(x,y) phi(x,y)^((n-3)/2)*dxPhi(x,y);
    dx = .001;
    L = max(appCtx.mesh.nodes(1,:));
    x = 0:dx:L;
    funcVal = zeros(size(x)); y_0 = 1;
    Linf = 0; L1 = 0; L2 = 0;
    
    for j=1:length(x)
        Linf = max(fun(x(j),y_0),Linf);
        L1 = L1 + abs(fun(x(j),y_0));
        L2 = L2 + fun(x(j),y_0)^2;
        funcVal(j) = fun(x(j),y_0);
    end
    L1 = L1*dx;
    L2 = sqrt(L2*dx);
    
    % get values of computed solution at quadrature points
    [xPoints, yPoints, uVal] = projectDOF(u, appCtx);
    
    xGrid = unique(sort(xPoints));
    yGrid = unique(sort(yPoints));
    [X,Y] = meshgrid(xGrid,yGrid);
    
    % plot fluid pressure
    funPlot = griddata(xPoints,yPoints,uVal.field{fFluidPres}(1,:)',X,Y);
    funPlot = funPlot(1,:);
    subplot(1,3,1), hold on
    h(i) = plot(xGrid,funPlot,'color',cc(i,:),'LineWidth', 3);
    view(90,90)
    xlabel('depth','FontSize', 15)
    title(' fluid pressure','FontSize', 17)
    
    % plot relative velocity
    funPlot = griddata(xPoints,yPoints,uVal.field{fRelVel}(1,:)',X,Y);
    funPlot = funPlot(1,:);
    subplot(1,3,2), hold on
    plot(xGrid,funPlot,'color',cc(i,:),'LineWidth', 3);
    view(90,90)
    title(' relative velocity','FontSize', 17)
    
    % plot porosity
    phiVal = zeros(1,length(xGrid));
    for j = 1:length(xGrid)
        phiVal(j) = phi(xGrid(j),yGrid(1));
    end
    subplot(1,3,3), hold on
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
    subplot(1,3,1), hold on
    funPlot = griddata(xPoints,yPoints,uExactVal.field{fFluidPres}(1,:)',X,Y);
    funPlot = funPlot(1,:);
    h(lenPhiVec+1) = plot(xGrid, funPlot,'color',cc(lenPhiVec+1,:),'LineWidth', 3,'LineStyle','--');
    legendStr{lenPhiVec+1} = 'exact solution';
    
    % plot exact relative velocity
    subplot(1,3,2), hold on
    funPlot = griddata(xPoints,yPoints,uExactVal.field{fRelVel}(1,:)',X,Y);
    funPlot = funPlot(1,:);
    plot(xGrid, funPlot,'color',cc(lenPhiVec+1,:),'LineWidth', 3,'LineStyle','--');
end

legend(h, legendStr);
