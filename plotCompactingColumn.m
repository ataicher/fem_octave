function outputInfo = plotCompactingColumn(u, appCtx)

FORMULATION = appCtx.FORMULATION;
Theta = appCtx.Theta;
phi = appCtx.phi;
EXISTEXACTSOL = appCtx.EXISTEXACTSOL;

% get values of computed and exact (if exists) solution at quadrature points
[xPoints, yPoints, uVal] = projectDOF(u, appCtx);
if EXISTEXACTSOL
    [xPoints, yPoints, uExactVal] = evalExactSolution(appCtx);
end

figure, clf,

xGrid = unique(sort(xPoints));
yGrid = unique(sort(yPoints));
[X,Y] = meshgrid(xGrid,yGrid);
outputInfo.xGrid = xGrid';

% plot porosity
phiVal = zeros(1,length(xGrid));
for i = 1:length(xGrid)
    phiVal(i) = phi(xGrid(i),yGrid(1));
end

outputInfo.phi = phiVal;

subplot(1,5,1)
plot(xGrid, phiVal, 'LineWidth', 4)
view(90,90)
ax = gca;
set(ax,'Box','on','FontSize',15);
title('porosity','FontSize', 15)
xlabel('depth','FontSize', 20)

% get field numbers
DarcyVelField = 1;
relVelField = 1;
fPresField = 2;
sfPresField = 2;
matVelField = 3;
mPresField = 4;


% plot Darcy velocity and matrix velocity
DarcyVel = griddata(xPoints,yPoints,uVal.field{DarcyVelField}(1,:)',X,Y);
if FORMULATION == 3
    DarcyVel = DarcyVel(1,:).*(phiVal.^(1+Theta));
elseif FORMULATION == 4
    relVelDOF = u(appCtx.field(DarcyVelField).startDOF:appCtx.field(DarcyVelField).endDOF);
    DarcyVelDOF = relVelDOF.*(appCtx.degenerate.phi.e.^(1+Theta));
    uTmp = u;
    uTmp(appCtx.field(DarcyVelField).startDOF:appCtx.field(DarcyVelField).endDOF) = DarcyVelDOF;
    [xPoints, yPoints, uValTmp] = projectDOF(uTmp, appCtx);
    DarcyVel = griddata(xPoints,yPoints,uValTmp.field{DarcyVelField}(1,:)',X,Y);
    DarcyVel = DarcyVel(1,:);
else
    DarcyVel = DarcyVel(1,:);
end
matVel = griddata(xPoints,yPoints,uVal.field{matVelField}(1,:)',X,Y);
matVel = matVel(1,:);

if EXISTEXACTSOL
    DarcyVelE = griddata(xPoints,yPoints,uExactVal.field{DarcyVelField}(1,:)',X,Y);
    if FORMULATION == 3 || FORMULATION == 4
        DarcyVelE = DarcyVelE(1,:).*(phiVal.^(1+Theta));
    else
        DarcyVelE = DarcyVelE(1,:);
    end
    matVelE = griddata(xPoints,yPoints,uExactVal.field{matVelField}(1,:)',X,Y);
    matVelE = matVelE(1,:);
end

outputInfo.u = DarcyVel;
outputInfo.v_m = matVel;
if EXISTEXACTSOL
    outputInfo.uE = DarcyVelE;
    outputInfo.v_mE = matVelE;
end

subplot(1,5,2); hold on
if EXISTEXACTSOL
    plot(xGrid, DarcyVelE, 'LineWidth', 3,'color','red')
    plot(xGrid, matVelE, 'LineWidth', 3,'color','magenta')
end
plot(xGrid, DarcyVel, 'LineWidth', 3,'color','black','LineStyle','--')
plot(xGrid, matVel, 'LineWidth', 3,'color','blue','LineStyle','--')
view(90,90)
ax = gca;
set(ax,'Box','on','FontSize',15);
set(ax,'XTick',[])
title('Darcy and matrix velocity','FontSize', 15)
if EXISTEXACTSOL
    legend('uE','v_mE','u','v_m','Location','NorthEast')
else
    legend('u','v_m','Location','NorthEast')
end

% plot fluid and matrix pressure potentials
if FORMULATION == 4
    fPresDOF = u(appCtx.field(fPresField).startDOF:appCtx.field(fPresField).endDOF);
    fPresDOF = fPresDOF./sqrt(appCtx.degenerate.phi.c);
    uTmp = u;
    uTmp(appCtx.field(fPresField).startDOF:appCtx.field(fPresField).endDOF) = fPresDOF;
    [xPoints, yPoints, uValTmp] = projectDOF(uTmp, appCtx);
    fPres = griddata(xPoints,yPoints,uValTmp.field{fPresField}(1,:)',X,Y);
    fPres(fPres==Inf) = 0;
    fPres(isnan(fPres)) = 0;
else
    fPres = griddata(xPoints,yPoints,uVal.field{fPresField}(1,:)',X,Y);
    fPres = fPres(1,:);
end
mPres = griddata(xPoints,yPoints,uVal.field{mPresField}(1,:)',X,Y);
mPres = mPres(1,:);

if EXISTEXACTSOL
    fPresE = griddata(xPoints,yPoints,uExactVal.field{fPresField}(1,:)',X,Y);
    if FORMULATION == 4
        fPresE = fPresE(1,:)./sqrt(phiVal);
        fPresE(fPresE==Inf) = 0;
        fPresE(isnan(fPresE)) = 0;
    else
        fPresE = fPresE(1,:);
    end
    mPresE = griddata(xPoints,yPoints,uExactVal.field{mPresField}(1,:)',X,Y);
    mPresE = mPresE(1,:);
    
end

outputInfo.q_f = fPres;
outputInfo.q_m = mPres;
if EXISTEXACTSOL
    outputInfo.q_fE = fPresE;
    outputInfo.q_mE = mPresE;
end

subplot(1,5,3); hold on
if EXISTEXACTSOL
    plot(xGrid, fPresE, 'LineWidth', 3,'color','red')
    plot(xGrid, mPresE, 'LineWidth', 3,'color','magenta')
end
plot(xGrid, fPres, 'LineWidth', 3,'color','black','LineStyle','--')
plot(xGrid, mPres, 'LineWidth', 3,'color','blue','LineStyle','--')
view(90,90)
ax = gca;
set(ax,'Box','on','FontSize',15);
set(ax,'XTick',[])
title('fluid and matrix pres. pot.','FontSize', 15)

if EXISTEXACTSOL
    legend('q_fE','q_mE','q_f','q_m','Location','NorthEast')
else
    legend('q_f','q_m','Location','NorthEast')
end

% plot relative velocity
if FORMULATION == 3 || FORMULATION == 4
    relVel = griddata(xPoints,yPoints,uVal.field{relVelField}(1,:)',X,Y);
    relVel = relVel(1,:);
    
    if EXISTEXACTSOL
        relVelE = griddata(xPoints,yPoints,uExactVal.field{relVelField}(1,:)',X,Y);
        relVelE = relVelE(1,:);
    end
    
    outputInfo.v_r = relVel;
    if EXISTEXACTSOL
        outputInfo.v_rE = relVelE;
    end
    
    subplot(1,5,4); hold on
    if EXISTEXACTSOL
        plot(xGrid, relVelE, 'LineWidth', 3,'color','red')
    end
    plot(xGrid, relVel, 'LineWidth', 3,'color','black','LineStyle','--')
    ax = gca;
    set(ax,'Box','on','FontSize',15);
    view(90,90)
    set(ax,'XTick',[])
    title('relative velocity','FontSize', 15)
    
    if EXISTEXACTSOL
        legend('v_rE','v_r','Location','NorthEast')
    else
        legend('v_r','Location','NorthEast')
    end
end

% plot scaled pressure potential
if FORMULATION == 4
    sfPres = griddata(xPoints,yPoints,uVal.field{sfPresField}(1,:)',X,Y);
    sfPres = sfPres(1,:);
    
    if EXISTEXACTSOL
        sfPresE = griddata(xPoints,yPoints,uExactVal.field{sfPresField}(1,:)',X,Y);
        sfPresE = sfPresE(1,:);
    end
    
    outputInfo.sq_f = sfPres;
    if EXISTEXACTSOL
        outputInfo.sq_fE = sfPresE;
    end
    
    subplot(1,5,5); hold on
    if EXISTEXACTSOL
        plot(xGrid, sfPresE, 'LineWidth', 3,'color','red')
    end
    plot(xGrid, sfPres, 'LineWidth', 3,'color','black','LineStyle','--')
    ax = gca;
    set(ax,'Box','on','FontSize',15);
    view(90,90)
    set(ax,'XTick',[])
    title('scaled fluid pres. pot.','FontSize', 15)
    if EXISTEXACTSOL
        legend('~q_fE','~q_f','Location','NorthEast')
    else
        legend('~q_f','Location','NorthEast')
    end
end