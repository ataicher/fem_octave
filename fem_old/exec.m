function [u, norms, errors, appCtx] = exec(input, varargin)

%% SET PATH TO LOCAL MATLAB DIRECTORY AND SET FIGURES TO SINGLE WINDOW-----
restoredefaultpath;
addpath(genpath(pwd));
set(0,'DefaultFigureWindowStyle','docked')

%% PROCESS INPUT FILE
% get optional command line input
% 4. nx
% 5. ny
% 6. set porosity
% 7. set minimum porosity
% 8. pick formulation
if nargin > 4
    nx = varargin{3}; ny = varargin{4};
    if nargin > 5
        SETPHI = varargin{5};
        if nargin > 6
            phi_m = varargin{6};
            if nargin > 7
                phi_p = varargin{7};
                if nargin > 8
                    FORMULATION = varargin{8};
                end
            end
        end
    end
end
input();

%% COLLECT PROBLEM INFORMATION---------------------------------------------
clear appCtx;
appCtx.DEBUGLEVEL = 0;
appCtx.OUTPUTMODE = 0;
if nargin > 1
    appCtx.DEBUGLEVEL = varargin{1};
    if nargin > 2
        appCtx.OUTPUTMODE = varargin{2};
    end
end
assembleAppCtx;
numFields = appCtx.numFields;
h = appCtx.h;
relTol = appCtx.NewtonRelTol;
absTol = appCtx.NewtonAbsTol;
maxIter = appCtx.NewtonMaxIter;
EXISTEXACTSOL = appCtx.EXISTEXACTSOL;
DEBUGLEVEL= appCtx.DEBUGLEVEL;
OUTPUTMODE = appCtx.OUTPUTMODE;

%% OUTPUT PROBLEM INFO-----------------------------------------------------
if DEBUGLEVEL > 1
    figure, clf, plotMesh(appCtx);
    figure, clf, plotBndryCondition(appCtx);
    outputProblemInfo(appCtx);
    if DEBUGLEVEL > 2
        plotBasisFuncs(appCtx);
    end
end

%% INITIAL GUESS-----------------------------------------------------------
% initial guess determined by essential condition and nullspace
if DEBUGLEVEL > 0
    fprintf('\n...calculating initial guess')
end
u = lift(appCtx);
if DEBUGLEVEL > 0
    fprintf('...done\n')
end

% display initial guess
if DEBUGLEVEL > 1
    [xPoints, yPoints, uVal] = projectDOF(u, appCtx);
    figure, clf, plotFun(xPoints, yPoints, uVal, appCtx);
    suptitle('Initial guess')
end

%% INITIAL RESIDUAL--------------------------------------------------------
% compute residual
if DEBUGLEVEL > 0
    fprintf('\n...calculating initial residual')
end
r = computeResidual(u,appCtx);
if DEBUGLEVEL > 0
    fprintf('...done\n')
end

% initial residual norm
residualErrorInitial = residualNorm(r, appCtx);
residualError = residualErrorInitial;
if DEBUGLEVEL > 0
    fprintf('\nInitial residual error L2 Norm:\n');
    for f=1:numFields
        name = appCtx.field(f).name;
        fprintf('   %-30.30s: %.4e\n', name, residualError(f));
    end
    fprintf('\n')
end

%% NEWTON METHOD-----------------------------------------------------------
if DEBUGLEVEL > 1
    figure,
end
iter = 0;
while (iter == 0) || (norm(residualError) > relTol*norm(residualErrorInitial) ...
        && norm(residualError) > absTol && iter <= maxIter)
    
    if DEBUGLEVEL > 0
        fprintf('\nNewton iteration: %d\n', iter)
    end
    
    % check finite difference Jacobian converges to provided Jacobian form
    if DEBUGLEVEL > 2 && h == 0
        fprintf('\n   ...checking Jacobian form')
        figure, clf, checkJacobianForm(appCtx);
        fprintf('...done\n')
    end
    
    % compute Jacobian matrix using user defined Jacobian form or finite
    % difference scheme
    if DEBUGLEVEL > 0
        fprintf('\n   ...calculating Jacobian')
    end
    Jac = computeJacobian(u,appCtx);
    if DEBUGLEVEL > 0
        fprintf('...done\n')
    end
    
    % compute Jacobian condition number
    if DEBUGLEVEL > 0
        fprintf('\n   ...calculating Jacobian condition number')
    end
    if issparse(Jac)
        appCtx.condJac = condest(Jac);
    else
        appCtx.condJac = 1/rcond(Jac);
    end
    if DEBUGLEVEL > 0
        fprintf('...done\n')
        fprintf('\n   Estimated Jacobian condition number: %e\n',appCtx.condJac);
    end
    
    % compute next iterate
    if DEBUGLEVEL > 0
        fprintf('\n   ...solving linear system');
    end
    u = u - Jac\r;
    if DEBUGLEVEL > 0
        fprintf('...done\n');
    end
    
    % view current iterate evaluated at quadrature points
    if DEBUGLEVEL > 1
        [xPoints, yPoints, uVal] = projectDOF(u, appCtx);
        clf, plotFun(xPoints, yPoints, uVal,appCtx)
        suptitle('current iterate of Newton method')
        pause(2)
    end
    
    % comupte new residual
    if DEBUGLEVEL > 0
        fprintf('\n   ...updating residual')
    end
    r = computeResidual(u,appCtx);
    if DEBUGLEVEL > 0
        fprintf('...done\n')
    end
    
    % display residual error
    residualError = residualNorm(r, appCtx);
    if DEBUGLEVEL > 0
        fprintf('\n   Residual error L2 Norm:\n');
        for f=1:numFields
            name = appCtx.field(f).name;
            fprintf('      %-30.30s: %.4e\n', name, residualError(f));
        end
        fprintf('\n')
    end
    
    % update Newton iteration
    iter = iter + 1;
end

% check if max number of iteration has been reached
if norm(residualError) > relTol*norm(residualErrorInitial) ...
        || norm(residualError) > absTol
    fprintf('\nNewton Method failed to converge: max number of Newton iteration reached\n')
elseif DEBUGLEVEL > 0
    fprintf('\nNewton Method converged in %d iteration(s)\n', iter)
end

%% COMPUTED SOLUTION-------------------------------------------------------
if DEBUGLEVEL > 1
    suptitle('Computed solution evaluated on quad points')
elseif DEBUGLEVEL == 1
    [xPoints, yPoints, uVal] = projectDOF(u, appCtx);
    figure, clf, plotFun(xPoints, yPoints, uVal,appCtx)
    suptitle('Computed solution evaluated on quad points')
end
    
%% JACOBIAN SPYPLOT--------------------------------------------------------
if DEBUGLEVEL > 1
    figure, clf, plotJacobianSpyPlot(Jac,appCtx);
end

%% STREAM FUNCTION---------------------------------------------------------
if DEBUGLEVEL > 1
    if appCtx.EXISTSTREAMFUN
        [xPoints, yPoints, streamFun] = computeStreamFunction(u,appCtx);
        figure, clf, plotStreamLines(xPoints, yPoints, appCtx, streamFun);
    else
        figure, clf, plotStreamLines(xPoints, yPoints, appCtx, u);
    end
    suptitle('Stream lines')
end

%% ERROR AND NORMS
if DEBUGLEVEL > 0
    fprintf('\n...calculating norms and errors')
end
if EXISTEXACTSOL
    [norms, errors] = calcNormAndRelError(u,appCtx);
else
    norms = calcNormAndRelError(u,appCtx);
end
if DEBUGLEVEL > 0
    fprintf('...done\n')
end

if EXISTEXACTSOL
    
    if DEBUGLEVEL > 0
        % display exact solution
        [xPoints, yPoints, uExactVal] = evalExactSolution(appCtx);
        figure, clf, plotFun(xPoints, yPoints, uExactVal, appCtx)
        suptitle('Exact solution evaluated on quad points')
    end
    if DEBUGLEVEL > 1
        % display differnce between exact and computed solution
        for f=1:numFields
            errorVal.field{f} = abs(uExactVal.field{f} - uVal.field{f});
        end
        figure, clf, plotFun(xPoints, yPoints, errorVal, appCtx)
        suptitle('Error evaluated on quad points')
    end
    
            uExactDOF = projectExactSolution(appCtx);
            
    if DEBUGLEVEL > 1
        % display interploant of exact solution
        [xPoints, yPoints, uExactVal] = projectDOF(uExactDOF, appCtx);
                figure, clf, plotFun(xPoints, yPoints, uExactVal, appCtx)
        suptitle('Interpolant of exact solution evaluated on quad points')
    end
    if DEBUGLEVEL > 0
        % display interpolation error between exact and computed solution
        [xPoints, yPoints, uExactVal] = projectDOF(uExactDOF, appCtx);
        for f=1:numFields
            errorVal.field{f} = abs(uExactVal.field{f} - uVal.field{f});
        end
        figure, clf, plotFun(xPoints, yPoints, errorVal, appCtx)
        suptitle('Interpolation error evaluated on quad points')
    end
    
end

%% RESIDUAL ERROR OF COMPUTED SOLUTION
for f=1:numFields
    errors(f).residualError = residualError(f);
end
if DEBUGLEVEL > 1
    [xPoints, yPoints, rVal] = projectDOF(r, appCtx);
    figure, clf, plotFun(xPoints, yPoints, rVal, appCtx)
    suptitle('Residual of computed solution evaluated at quadrature points')
end

% output errors and norms
if OUTPUTMODE
    outputNormsErrors(norms, errors, appCtx);
end
