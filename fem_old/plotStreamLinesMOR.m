function plotStreamLinesMOR(u, appCtx)
figure,clf
[xPoints, yPoints, streamFun] = computeStreamFunction(u,appCtx);

xGrid = unique(sort(xPoints));
yGrid = unique(sort(yPoints));
xMin = min(xGrid); xMax = max(xGrid);
yMin = min(yGrid); yMax = max(yGrid);
[X,Y] = meshgrid(xGrid,yGrid);

for f=1:4
    name = appCtx.field(f).name;
    if strcmp(name,'relative velocity') || strcmp(name, 'Darcy velocity')
        relVelocityField = f;
    elseif strcmp(name, 'matrix velocity')
        matVelocityField = f;
    end
end
streamFun{relVelocityField} = streamFun{relVelocityField} + streamFun{matVelocityField};
StreamFun = griddata(xPoints,yPoints,streamFun{relVelocityField},X,Y);
v = StreamFun(:,end);
contour(X,Y,StreamFun,v,'--','LineWidth',2)
hold on
StreamFun = griddata(xPoints,yPoints,streamFun{matVelocityField},X,Y);
v = StreamFun(:,end);
contour(X,Y,StreamFun,v,'LineWidth',2)
xlabel('depth','FontSize',20), ylabel('x','FontSize',20)
axis([xMin xMax yMin yMax])
title('Stream function. Matrix - solid lines. Fluid - dashed lines')

