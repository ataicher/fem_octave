function outputNormsErrors(norms,errors,appCtx)

DEBUGLEVEL = appCtx.DEBUGLEVEL;
EXISTEXACTSOL = appCtx.EXISTEXACTSOL;
EXISTEXACTSOLGRAD = appCtx.EXISTEXACTSOLGRAD;
numFields = appCtx.numFields;

if DEBUGLEVEL > 0
    fprintf('\n------------------------------------------------------------\n')
end

for f=1:numFields
    name = appCtx.field(f).name;
    fprintf('\n\n%s:\n', name);
    if appCtx.field(f).CALCHDIVNORM
        if EXISTEXACTSOL
            if EXISTEXACTSOLGRAD
                fprintf('%14.14s | %14.14s | %14.14s | %14.14s | %14.14s | %14.14s |\n', 'L2 norm', 'L2 error', 'L2 int-error', 'H(div) norm', 'H(div) error', 'H(div) int-error');
                fprintf('-----------------------------------------------------------------------------------------------------\n');
                fprintf('%14.4e | %14.4e | %14.4e | %14.4e | %14.4e | %14.4e |\n', norms(f).L2, errors(f).L2, errors(f).L2Int, norms(f).HDiv, errors(f).HDiv, errors(f).HDivInt);
            else
                fprintf('%14.14s | %14.14s | %14.14s | %14.14s |\n', 'L2 norm', 'L2 error', 'L2 int-error', 'HDiv norm');
                fprintf('--------------------------------------------------------\n');
                fprintf('%14.4e | %14.4e | %14.4e | %14.4e |\n', norms(f).L2, errors(f).L2, errors(f).L2Int, norms(f).HDiv);
            end
        else
            fprintf('%14.14s | %14.14s |\n', 'L2 norm', 'HDiv norm');
            fprintf('---------------------------------\n');
            fprintf('%14.4e | %14.4e |\n', norms(f).L2, norms(f).HDiv);
        end
        
    elseif appCtx.field(f).CALCH1NORM
        if EXISTEXACTSOL
            if EXISTEXACTSOLGRAD
                fprintf('%14.14s | %14.14s | %14.14s | %14.14s | %14.14s | %14.14s |\n', 'L2 norm', 'L2 error', 'L2 int-error', 'H1 norm', 'H1 error', 'H1 int-error');
                fprintf('-----------------------------------------------------------------------------------------------------\n');
                fprintf('%14.4e | %14.4e | %14.4e | %14.4e | %14.4e | %14.4e |\n', norms(f).L2, errors(f).L2, errors(f).L2Int, norms(f).H1, errors(f).H1, errors(f).H1Int);
            else
                fprintf('%14.14s | %14.14s | %14.14s | %14.14 |s\n', 'L2 norm', 'L2 error', 'L2 int-error', 'H1 norm');
                fprintf('-----------------------------------------------------\n');
                fprintf('%14.4e | %14.4e | %14.4e | %14.4e |\n', norms(f).L2, errors(f).L2, errors(f).L2Int, norms(f).H1);
            end
        else
            fprintf('%14.14s | %14.14s |\n', 'L2 norm', 'H1 norm');
            fprintf('---------------------------------\n');
            fprintf('%14.4e | %14.4e |\n', norms(f).L2, norms(f).H1);
        end
    else
        if EXISTEXACTSOL
            fprintf('%14.14s | %14.14s | %14.14s |\n', 'L2 norm', 'L2 error', 'L2 int-error');
            fprintf('--------------------------------------------------\n');
            fprintf('%14.4e | %14.4e | %14.4e |\n', norms(f).L2, errors(f).L2, errors(f).L2Int);
        else
            fprintf('%14.14s |\n', 'L2 norm');
            fprintf('----------------\n');
            fprintf('%14.4e |\n', norms(f).L2);
        end
    end
end
fprintf('\n');