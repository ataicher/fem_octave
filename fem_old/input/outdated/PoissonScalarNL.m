% POISSONSCALARNL
% Solve the PDE
%   -div((1+mu*u^2)*gradU) = f
% for all v in V
%    ((1+mu*u^2)*gradU , gradV) - (f , v) - <(1+mu*u^2)*gradU.*n , v> = 0

%% MESH--------------------------------------------------------------------
% uniform quad mesh on a rectangular domain
[nodes, edges, cells] = uniformQuadMesh(10,10,[-1 1],[-1 1]);
% [nodes, edges, cells] = QuadMesh([-1 -.9 -.8 -.4  1],[-1 .4 .8 .9 1]);

%% VARIABLES---------------------------------------------------------------
name.field{1} = 'pressure';
numComp.field{1} = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 2;

%% FINITE ELEMENT BASIS----------------------------------------------------
basisType.field{1} = @P1Quad;

%% MANUFACTURED SOLUTION---------------------------------------------------
uE.field{1}{1} = @(x,y) x*y;

[gradUE, grad2UE] = symDerivatives(uE);
U = uE.field{1}{1}; 
dxU = gradUE.field{1}{1,1}; dyU = gradUE.field{1}{1,2};
dxxU = grad2UE.field{1}{1,1,1}; dyyU = grad2UE.field{1}{1,2,2};

%% VARIATIONAL FORM--------------------------------------------------------
% nonlinearity
mu = 1;
% f = -div((1+mu*u^2)*gradU) = -2*mu*u*(dxu^2 + dyu^2) - (1 + mu*u^2)*(dxxu + dyyu)
f = @(x,y) -2*mu*U(x,y)*(dxU(x,y)^2 + dyU(x,y)^2) - ...
    (1 + mu*U(x,y)^2)*(dxxU(x,y) + dyyU(x,y));
% 
varForm.field{1}.v{1} = @(u, gradU, x, y) -f(x,y);
% 
varForm.field{1}.gradV{1,1} = @(u, gradU, x, y) (1 + mu*u.field{1}^2)*gradU.field{1}(1);
varForm.field{1}.gradV{1,2} = @(u, gradU, x, y) (1 + mu*u.field{1}^2)*gradU.field{1}(2);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0;
% ((1+mu*u^2)*gradW , gradV) 
JacForm.field{1,1}.gradW{1,1}.gradV{1,1} = @(u,gradU,x,y) 1 + mu*u.field{1}^2;
JacForm.field{1,1}.gradW{1,2}.gradV{1,2} = @(u,gradU,x,y) 1 + mu*u.field{1}^2;
% (2*mu*u*gradU*w , gradV)
JacForm.field{1,1}.w{1}.gradV{1,1} = @(u,gradU,x,y) 2*mu*u.field{1}*gradU.field{1}(1);
JacForm.field{1,1}.w{1}.gradV{1,2} = @(u,gradU,x,y) 2*mu*u.field{1}*gradU.field{1}(2);

%% BOUNDARY CONDITIONS-----------------------------------------------------
% bottom edge
loc = @(x,y) y+1;
n = [0 -1];
alpha = 0;
beta = 1;
eta = @(x,y) -alpha*(1 + mu*U(x,y)^2)*(dxU(x,y)*n(1) + dyU(x,y)*n(2)) + ...
    beta*U(x,y);
BC.field{1}.bndry{1} = assembleBoundary(loc,alpha,beta,eta);

% right edge
loc = @(x,y) x-1;
n = [1 0];
alpha = 0;
beta = 1;
eta = @(x,y) -alpha*(1 + mu*U(x,y)^2)*(dxU(x,y)*n(1) + dyU(x,y)*n(2)) + ...
    beta*U(x,y);
BC.field{1}.bndry{2} = assembleBoundary(loc,alpha,beta,eta);

% top edge
loc = @(x,y) y-1;
n = [0 1];
alpha = 0;
beta = 1;
eta = @(x,y) -alpha*(1 + mu*U(x,y)^2)*(dxU(x,y)*n(1) + dyU(x,y)*n(2)) + ...
    beta*U(x,y);
BC.field{1}.bndry{3} = assembleBoundary(loc,alpha,beta,eta);

% left edge
loc = @(x,y) x+1;
n = [-1 0];
alpha = 0;
beta = 1;
eta = @(x,y) -alpha*(1 + mu*U(x,y)^2)*(dxU(x,y)*n(1) + dyU(x,y)*n(2)) + ...
    beta*U(x,y);
BC.field{1}.bndry{4} = assembleBoundary(loc,alpha,beta,eta);

%% SOLVER------------------------------------------------------------------
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 10;

%% NULLSPACE---------------------------------------------------------------
% pick fields that have a constant nullSpace
% nullSpace.field{1} = 'true';

%% DISPLAY AND DEBUG OPTIONS-----------------------------------------------
displayMesh = 'true';
displayInitialGuess = 'true';
displayResidualError = 'true';
% displayCurrentNewtonIterate = 'true';
displayComputedSolution = 'true';
exactSolutionExist = 'true';
displayExactSolution = 'true';
% displayStreamLines = 'true';
% displayJacobianSpyPlot = 'true';