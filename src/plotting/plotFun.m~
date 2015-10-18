function plotFun(xPoints, yPoints, funcVal, appCtx)

numFields = appCtx.numFields;
maxNumComp = appCtx.maxNumComp;

xGrid = unique(sort(xPoints));
yGrid = unique(sort(yPoints));
[X,Y] = meshgrid(xGrid,yGrid);

hTabGroup = uitabgroup; drawnow;
warning('off', 'MATLAB:uitabgroup:OldVersion');
for f=1:numFields
    numComp = appCtx.field(f).numComp;
    name = appCtx.field(f).name;
    
    subPlotInd = maxNumComp*(f-1);
    for comp=1:numComp
        funPlot = griddata(xPoints,yPoints,funcVal.field{f}(comp,:)',X,Y);
        
        % ind = find(Y(:,1)<-1);
        % Y = Y(ind,:);
        % X = X(ind,:);
        % funPlot = funPlot(ind,:);
        subPlotInd = subPlotInd + 1;
        
        if numComp > 1
            titleName = [name, ' comp ',num2str(comp)];
        else
            titleName = name;
        end
        
        subPlotInd = subPlotInd + 1;
        tab(subPlotInd) = uitab(hTabGroup, 'title',titleName);
        a = axes('parent', tab(subPlotInd));
        surf(X,Y,funPlot)
        caxis([-10 10])
        colorbar
        xlabel('x')
        ylabel('y')
        title(titleName,'FontSize',15)
    end
end
% uicontrol(tab(end), 'String','Close', 'Callback','close(gcbf)');