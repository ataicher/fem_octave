function convRate = calcConvRates(input, nRef, nxStart, nyStart)

%% SET PATH TO LOCAL MATLAB DIRECTORY AND SET FIGURES TO SINGLE WINDOW-----
restoredefaultpath;
addpath(genpath(pwd));
set(0,'DefaultFigureWindowStyle','docked')

% get problem info
input();
numFields = length(field);

% no refinement if in x-dir if nxStart = 0, y-dir if nystart = 0
if nxStart
    nx = nxStart; 
else
    nx = 1;
end
if nyStart
    ny = nyStart;
else
    ny = 1;
end
nnx = [];
nny = [];
for ref=1:nRef
    
    % display nx and ny
    fprintf('\nsolving for nx = %d, ny = %d\n',nx,ny);
    
    nnx = [nnx nx];
    nny = [nny ny];
    
    % run driver for this refinement
    [~, ~, errors, appCtx] = exec(input,0,0,nx,ny);
    
    for f = 1:numFields
        name = appCtx.field(f).name;
        
        % get errors and check that they are nonzero
        if abs(errors(f).L2)>1e-10
            RefErrors(f,ref).L2 = errors(f).L2;
            if ref > 1
                convRate(f,ref).L2 = log(RefErrors(f,ref-1).L2/RefErrors(f,ref).L2)/log(2);
            else
                convRate(f,ref).L2 = 0;
            end
        else
            fprintf('L2 error converged within tolerance 1e-10 for field %s\n', name)
        end
        
        if abs(errors(f).L2Int)>1e-10
            RefErrors(f,ref).L2Int = errors(f).L2Int;
            if ref > 1 
                convRate(f,ref).L2Int = log(RefErrors(f,ref-1).L2Int/RefErrors(f,ref).L2Int)/log(2);
            else
                convRate(f,ref).L2Int = 0;
            end
        else
            fprintf('L2 interpolation error converged within tolerance 1e-10 for %s\n', name)         
        end
        
        if appCtx.EXISTEXACTSOLGRAD
            if appCtx.field(f).CALCHDIVNORM
                if abs(errors(f).HDiv)>1e-10
                    RefErrors(f,ref).HDiv = errors(f).HDiv;
                    if ref > 1
                        convRate(f,ref).HDiv = log(RefErrors(f,ref-1).HDiv/RefErrors(f,ref).HDiv)/log(2);
                    else
                        convRate(f,ref).HDiv = 0;
                    end
                else
                    fprintf('H(div) error converged within tolerance 1e-10 for field %s\n', name)
                end
                
                if abs(errors(f).HDivInt)>1e-10
                    RefErrors(f,ref).HDivInt = errors(f).HDivInt;
                    if ref > 1
                        convRate(f,ref).HDivInt = log(RefErrors(f,ref-1).HDivInt/RefErrors(f,ref).HDivInt)/log(2);
                    else
                        convRate(f,ref).HDivInt = 0;
                    end
                else
                    fprintf('H(div) interpolation error converged within tolerance 1e-10 for field %s\n', name)
                end
            end
            
            if appCtx.field(f).CALCH1NORM
                if abs(errors(f).H1)>1e-10
                    RefErrors(f,ref).H1 = errors(f).H1;
                    if ref > 1
                        convRate(f,ref).H1 = log(RefErrors(f,ref-1).H1/RefErrors(f,ref).H1)/log(2);
                    else
                        convRate(f,ref).H1 = 0;
                    end
                else
                    fprintf('H1 error converged within tolerance 1e-10 for field %s\n', name)
                end
                
                if abs(errors(f).H1Int)>1e-10
                    RefErrors(f,ref).H1Int = errors(f).H1Int;
                    if ref > 1
                        convRate(f,ref).H1Int = log(RefErrors(f,ref-1).H1Int/RefErrors(f,ref).H1Int)/log(2);
                    else
                        convRate(f,ref).H1Int = 0;
                    end
                else
                    fprintf('H1 interpolation error converged within tolerance 1e-10 for field %s\n', name)
                end
            end
        end
        
    end
    
    % refine mesh
    if nxStart && nyStart
        nx=2*nx; ny = 2*ny;
    elseif nxStart
        nx = 2*nx;
    elseif nyStart
        ny = 2*ny;
    end
end

for f = 1:numFields
    name = appCtx.field(f).name;
    fprintf('\n%s\n', name);
    
    if appCtx.field(f).CALCH1NORM && appCtx.EXISTEXACTSOLGRAD
        
        fprintf('%4.12s | %4.12s | %13.12s | %8.12s | %13.12s | %8.12s | %13.12s | %8.12s | %13.12s | %8.12s\n', ...
            'nx', 'ny', 'L2 error', 'rate', 'H1 error', 'rate', 'L2 int-error', 'rate', 'H1 int-error', 'rate');
        fprintf('-----------------------------------------------------------------------------------------------------------------------\n');
        for ref=1:nRef
            fprintf('%4d | %4d | %13.4d | %8.4f | %13.4d | %8.4f | %13.4d | %8.4f | %13.4d | %8.4f\n', ...
                nnx(ref), nny(ref), RefErrors(f,ref).L2, convRate(f,ref).L2, RefErrors(f,ref).H1, convRate(f,ref).H1, RefErrors(f,ref).L2Int, convRate(f,ref).L2Int, RefErrors(f,ref).H1Int, convRate(f,ref).H1Int);
        end
        
    elseif appCtx.field(f).CALCHDIVNORM && appCtx.EXISTEXACTSOLGRAD
        
        fprintf('%4.12s | %4.12s | %13.12s | %8.12s | %13.12s | %8.12s | %13.12s | %8.12s | %13.12s | %8.12s\n', ...
            'nx', 'ny', 'L2 error', 'rate', 'H(div) error', 'rate', 'L2 int-error', 'rate', 'H(div) int-error', 'rate');
        fprintf('-----------------------------------------------------------------------------------------------------------------------\n');
        for ref=1:nRef
            fprintf('%4d | %4d | %13.4d | %8.4f | %13.4d | %8.4f | %13.4d | %8.4f | %13.4d | %8.4f\n', ...
                nnx(ref), nny(ref), RefErrors(f,ref).L2, convRate(f,ref).L2, RefErrors(f,ref).HDiv, convRate(f,ref).HDiv, RefErrors(f,ref).L2Int, convRate(f,ref).L2Int, RefErrors(f,ref).HDivInt, convRate(f,ref).HDivInt);
        end
        
    else
        
        fprintf('%4.12s | %4.12s | %13.12s | %8.12s | %13.12s | %8.12s\n', ...
            'nx', 'ny', 'L2 error', 'rate', 'L2 int-error', 'rate');
        fprintf('------------------------------------------------------------------\n');
        for ref=1:nRef
            fprintf('%4d | %4d | %13.4d | %8.4f | %13.4d | %8.4f\n', ...
                nnx(ref), nny(ref), RefErrors(f,ref).L2, convRate(f,ref).L2, RefErrors(f,ref).L2Int, convRate(f,ref).L2Int);
        end
    end
    
    fprintf('\n\n');
end
