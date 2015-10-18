% POISSONSCALAR
% Solve the PDE [-1,1]^2
%   -grad divU = f
% for all v in V
%    (gradU , gradV) - (f , v) - <gradU.*n , v> = 0

%% MESH--------------------------------------------------------------------
% uniform rectangular mesh on [-1 1]^2
[nodes, edges, cells] = uniformQuadMesh(1,1,[-1 1],[-1 1]);

%% VARIABLES---------------------------------------------------------------
name.field{1} = 'pressure';
numComp.field{1} = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 2;

%% FINITE ELEMENT BASIS----------------------------------------------------
% basisType.field{1} = @piecewiseConst;
basisType.field{1} = @P1Quad;

%% MANUFACTURED SOLUTION---------------------------------------------------
uE.field{1}{1} = @(x,y) 1+sin(3.14*y);
[gradUE, grad2UE] = symDerivatives(uE);
U = uE.field{1}{1}; 
dxU = gradUE.field{1}{1,1}; dyU = gradUE.field{1}{1,2};
dxxU = grad2UE.field{1}{1,1,1}; dyyU = grad2UE.field{1}{1,2,2};

%% VARIATIONAL FORM--------------------------------------------------------
f = @(x,y) -(dxxU(x,y) + dyyU(x,y));

% (gradU , gradV)
varForm.field{1}.gradV{1,1} = @(u, gradU, x, y) gradU.field{1}(1,1);
varForm.field{1}.gradV{1,2} = @(u, gradU, x, y) gradU.field{1}(1,2);

% - (f,v)
varForm.field{1}.v{1} = @(u, gradU, x, y) - f(x,y);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0;

% (gradU , gradV)
JacForm.field{1,1}.gradW{1,1}.gradV{1,1} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.gradW{1,2}.gradV{1,2} = @(u,gradU,x,y) 1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
% bottom edge
loc = @(x,y) y+1;
n = [0 -1];
alpha = 1;
beta = 0;
eta = @(x,y) -alpha*(dxU(x,y)*n(1) + dyU(x,y)*n(2)) + beta*U(x,y);
BC.field{1}.bndry{1} = assembleBoundary(loc,alpha,beta,eta);

% right edge
loc = @(x,y) x-1;
n = [1 0];
alpha = 1;
beta = 0;
eta = @(x,y) -alpha*(dxU(x,y)*n(1) + dyU(x,y)*n(2)) + beta*U(x,y);
BC.field{1}.bndry{2} = assembleBoundary(loc,alpha,beta,eta);

% top edge
loc = @(x,y) y-1;
n = [0 1];
alpha = 1;
beta = 0;
eta = @(x,y) -alpha*(dxU(x,y)*n(1) + dyU(x,y)*n(2)) + beta*U(x,y);
BC.field{1}.bndry{3} = assembleBoundary(loc,alpha,beta,eta);

% left edge
loc = @(x,y) x+1;
n = [-1 0];
alpha = 1;
beta = 0;
eta = @(x,y) -alpha*(dxU(x,y)*n(1) + dyU(x,y)*n(2)) + beta*U(x,y);
BC.field{1}.bndry{4} = assembleBoundary(loc,alpha,beta,eta);

%% SOLVER------------------------------------------------------------------
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 1;

%% NULLSPACE---------------------------------------------------------------
% pick fields that have a constant nullSpace
nullSpace.field{1} = 'true';

%% DISPLAY AND DEBUG OPTIONS-----------------------------------------------
% displayMesh = 'true';
% displayInitialGuess = 'true';
displayResidualError = 'true';
% plotResidualError = 'true';
% displayCurrentNewtonIterate = 'true';
displayComputedSolution = 'true';
exactSolutionExist = 'true';
displayExactSolution = 'true';
% displayStreamFunction = 'true';
% displayJacobianSpyPlot = 'true';