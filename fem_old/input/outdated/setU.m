% NONLINEAR POISSON - VECTOR FIELD
% Solve the PDE
%   u = f
%   u = u_0 on boundary
% for all v in V
%    (u , v) - (f , v) = 0
%


%% MESH--------------------------------------------------------------------
% uniform quad mesh on a rectangular domain
[nodes, edges, cells] = uniformQuadMesh(10,10,[-1 1],[-1 1]);
% [nodes, edges, cells] = QuadMesh([-1 -.9 -.8 -.4  1],[-1 .4 .8 .9 1]);

%% VARIABLES---------------------------------------------------------------
name.field{1} = 'u';
numComp.field{1} = 2;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 2;

%% FINITE ELEMENT BASIS----------------------------------------------------
basisType.field{1} = @P1Quad2;

%% MANUFACTURED SOLUTION---------------------------------------------------

%% VARIATIONAL FORM--------------------------------------------------------
uE.field{1}{1} = @(x,y) sin(y^2+x);
uE.field{1}{2} = @(x,y) x*y+cos(x*y);

[gradUE] = symDerivatives(uE);
Ux = uE.field{1}{1}; Uy = uE.field{1}{2};
dxUx = gradUE.field{1}{1,1}; dyUx = gradUE.field{1}{1,2}; dxUy = gradUE.field{1}{2,1}; dyUy = gradUE.field{1}{2,2};

%% VARIATIONAL FORM--------------------------------------------------------
fx = @(x,y) Ux(x,y);
fy = @(x,y) Uy(x,y);
% 
varForm.field{1}.v{1} = @(u, gradU, x, y) u.field{1}(1) - fx(x,y);
varForm.field{1}.v{2} = @(u, gradU, x, y) u.field{1}(2) - fy(x,y);

%% JACOBAIN FORM-----------------------------------------------------------
h = 0;
%
JacForm.field{1,1}.w{1}.v{1} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.w{2}.v{2} = @(u,gradU,x,y) 1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
alpha = [0 0];
beta = [1 1];

% bottom edge
loc = @(x,y) y+1;
n = [0 -1]; 
tau = [1 0];
eta{1} = @(x,y) beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2));
BC.field{1}.bndry{1} = assembleBoundary(loc,alpha,beta,eta);
% right edge
loc = @(x,y) x-1;
n = [1 0];
tau = [0 1];
eta{1} = @(x,y) beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2));
BC.field{1}.bndry{2} = assembleBoundary(loc,alpha,beta,eta);
% top edge
loc = @(x,y) y-1;
n = [0 1];
tau = [-1 0];
eta{1} = @(x,y) beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2));
BC.field{1}.bndry{3} = assembleBoundary(loc,alpha,beta,eta);
% left edge
loc = @(x,y) x+1;
n = [-1 0];
tau = [0 -1];
eta{1} = @(x,y) beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2));
BC.field{1}.bndry{4} = assembleBoundary(loc,alpha,beta,eta);

%% SOLVER
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 1;

%% NULLSPACE
% pick fields that have a constant nullSpace
nullSpace.field{2} = 'true';

%% DISPLAY AND DEBUG OPTIONS
% displayMesh = 'true';
% displayInitialGuess = 'true';
% displayResidualError = 'true';
% displayCurrentNewtonIterate = 'true';
displayComputedSolution = 'true';
exactSolutionExist = 'true';
displayExactSolution = 'true';
% displayStreamLines = 'true';
% displayJacobianSpyPlot = 'true';