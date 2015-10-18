function appCtx = assembleExactSolution(appCtx)

dim = appCtx.dim;
numFields = appCtx.numFields;

if isfield(appCtx.field, 'uExact')
    
    for f=1:numFields
        numComp = appCtx.field(f).numComp;
        
        if length(appCtx.field(f).uExact) < numComp
            error('myApp:argChk','exact solution for field %d undefined',f)
        end
        
        for comp = 1:numComp
            
            if isa(appCtx.field(f).uExact{comp},'function_handle') && nargin(appCtx.field(f).uExact{comp}) == 2
                %appCtx.field(f).uExact{comp} = expandFun(appCtx.field(f).uExact{comp});
            else
                error('myApp:argChk','exact solution for field %d undefined',f)
            end
            
        end
    end
    appCtx.EXISTEXACTSOL = 1;
    
    if isfield(appCtx.field, 'gradUExact')
        for f=1:numFields
            numComp = appCtx.field(f).numComp;
            
            if size(appCtx.field(f).gradUExact,1) < numComp ||  size(appCtx.field(f).gradUExact,2) < dim
                error('myApp:argChk','gradient of exact solution for field %d undefined',f)
            end
            
            for comp =1:numComp
                for d = 1:dim
                    if isa(appCtx.field(f).gradUExact{comp,d},'function_handle') && nargin(appCtx.field(f).gradUExact{comp,d}) == 2
                        %appCtx.field(f).gradUExact{comp,d} = expandFun(appCtx.field(f).gradUExact{comp,d});
                    else
                        error('myApp:argChk','gradient of exact solution for field %d undefined',f)
                    end
                end
            end
        end
        appCtx.EXISTEXACTSOLGRAD = 1;
    else
        appCtx.EXISTEXACTSOLGRAD = 0;
    end
    
else
    appCtx.EXISTEXACTSOL = 0;
    appCtx.EXISTEXACTSOLGRAD = 0;
end
