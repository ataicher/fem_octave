function appCtx = assembleNewton(relTol, absTol, maxIter, appCtx)

appCtx.NewtonRelTol = relTol;
appCtx.NewtonAbsTol = absTol;
appCtx.NewtonMaxIter = maxIter;