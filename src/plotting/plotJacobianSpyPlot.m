function plotJacobianSpyPlot(Jac, appCtx)

numFields = appCtx.numFields;
globalSize = appCtx.globalSize;

spy(abs(Jac) > 1e-10);
hold on
for f=1:numFields
    startDOF = appCtx.field(f).startDOF;
    endDOF = appCtx.field(f).endDOF;
    plot([startDOF startDOF],[0 globalSize],'r-')
    plot([0 globalSize],[startDOF startDOF],'r-')
    plot([endDOF endDOF],[0 globalSize],'r-')
    plot([0 globalSize],[endDOF endDOF],'r-')
end
suptitle('Jacobian SpyPlot')