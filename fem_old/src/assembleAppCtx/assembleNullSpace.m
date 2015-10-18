function appCtx = assembleNullSpace(appCtx)

numFields = appCtx.numFields;

for f=1:numFields
    if isempty(appCtx.field(f).nullSpace);
        appCtx.field(f).NULLSPACE = 0;
    else
        if strcmp(appCtx.field(f).nullSpace, 'true')
            nullDOF = appCtx.field(f).endDOF;
            appCtx.field(f).NULLSPACE = nullDOF;
        else
            appCtx.field(f).NULLSPACE = 0;
        end
    end
end

appCtx.field = rmfield(appCtx.field,'nullSpace');