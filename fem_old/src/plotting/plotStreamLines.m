function plotStreamLines(xPoints, yPoints, appCtx, varargin)

dim = appCtx.dim;
numFields = appCtx.numFields;
nodes = appCtx.mesh.nodes;
EXISTSTREAMFUN = appCtx.EXISTSTREAMFUN;

xGrid = unique(sort(xPoints));
yGrid = unique(sort(yPoints));
xMin = min(xGrid); xMax = max(xGrid);
yMin = min(yGrid); yMax = max(yGrid);
[X,Y] = meshgrid(xGrid,yGrid);

hTabGroup = uitabgroup; drawnow;
vectorFields = [];
for f=1:numFields
    numComp = appCtx.field(f).numComp;
    if numComp == dim
        %         numVectorFields = numVectorFields + 1;
        vectorFields = [vectorFields f];
    end
end
vectorFieldInd = 0;

if EXISTSTREAMFUN
    streamFun = varargin{1};
    for f=vectorFields
        name = appCtx.field(f).name;
        
        vectorFieldInd = vectorFieldInd + 1;
        StreamFun = griddata(xPoints,yPoints,streamFun{f},X,Y);
        
        %         ind = find(Y(:,1)<-1.4);
        %         Y = Y(ind,:);
        %         X = X(ind,:);
        %         StreamFun = StreamFun(ind,:);
        
        tab(vectorFieldInd) = uitab(hTabGroup, 'title',name);
        a = axes('parent', tab(vectorFieldInd));
        v = StreamFun(:,end);
        contour(X,Y,StreamFun,v)
        xlabel('x')
        ylabel('y')
        axis([xMin xMax yMin yMax])
        title(name,'FontSize',15)
    end
else
    u = varargin{1};
    numYCells = length(find(nodes(1,:) == max(nodes(1,:)))) - 1;
    [xPoints, yPoints, uVal] = projectDOF(u, appCtx);
    for f=vectorFields
        name = appCtx.field(f).name;
        
        vectorFieldInd = vectorFieldInd + 1;
        tab(vectorFieldInd) = uitab(hTabGroup, 'title',name);
        a = axes('parent', tab(vectorFieldInd));
        U = griddata(xPoints,yPoints,uVal.field{f}(1,:)',X,Y);
        V = griddata(xPoints,yPoints,uVal.field{f}(2,:)',X,Y);
        yStart = zeros(numYCells,1);
        for i=1:numYCells
            yStart(i) = (yMax-yMin)/(2*numYCells) + (i-1)*(yMax-yMin)/numYCells;
        end
        xStart = xMax*ones(size(yStart));
        streamline(X,Y,U,V,xStart,yStart);
        xlabel('x')
        ylabel('y')
        axis([xMin xMax yMin yMax])
        title(name,'FontSize',15)
        
    end
end
