% DARCY
%   Solve the PDE
%       u + grad p = f 
%            div u = g 
%   for all v in V
%       (u , v) - (p , divV) + <p , v.*n> - (f , v) = 0
%   for all q in Q
%       -(divU , q) + (g , q) = 0
%

%% MESH--------------------------------------------------------------------
% uniform quad mesh on a rectangular domain
[nodes, edges, cells] = uniformQuadMesh(20,20,[-1 1],[-1 1]);

%% VARIABLES---------------------------------------------------------------
name.field{1} = 'velocity';
numComp.field{1} = 2;
name.field{2} = 'pressure';
numComp.field{2} = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 2;

%% FINITE ELEMENT BASIS----------------------------------------------------
basisType.field{1} = @RT0Velocity;
basisType.field{2} = @piecewiseConst;

%% MANUFACTURED SOLUTION
uE.field{1}{1} = @(x,y) .5;
uE.field{1}{2} = @(x,y) 0;
uE.field{2}{1} = @(x,y) (1-x)/2;

[gradUE, grad2UE] = symDerivatives(uE);
Ux = uE.field{1}{1}; Uy = uE.field{1}{2};
dxUx = gradUE.field{1}{1,1}; dyUy = gradUE.field{1}{2,2};
P = uE.field{2}{1}; 
dxP = gradUE.field{2}{1,1}; dyP = gradUE.field{2}{1,2};

%% VARIATIONAL FORM--------------------------------------------------------
fx = @(x,y) -Ux(x,y) + dxP(x,y);
fy = @(x,y) -Uy(x,y) + dyP(x,y); 
g = @(x,y) dxUx(x,y) + dyUy(x,y);
% (u , v) - (f , v) 
varForm.field{1}.v{1} = @(u, gradU, x, y) u.field{1}(1) - fx(x,y);
varForm.field{1}.v{2} = @(u, gradU, x, y) u.field{1}(2) - fy(x,y);
% -(p , divV)
varForm.field{1}.gradV{1,1} = @(u, gradU, x, y) -u.field{2};
varForm.field{1}.gradV{2,2} = @(u, gradU, x, y) -u.field{2};
% -(divU, w) + (g , w)
varForm.field{2}.v{1} = @(u, gradU, x, y) -gradU.field{1}(1,1) - gradU.field{1}(2,2) + g(x,y);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0;
% -(w,v)
JacForm.field{1,1}.w{1}.v{1} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.w{2}.v{2} = @(u,gradU,x,y) 1;
% -(w,divV)
JacForm.field{2,1}.w{1}.gradV{1,1} = @(u,gradU,x,y) -1;
JacForm.field{2,1}.w{1}.gradV{2,2} = @(u,gradU,x,y) -1;
% -(divW,v)
JacForm.field{1,2}.gradW{1,1}.v{1} = @(u,gradU,x,y) -1;
JacForm.field{1,2}.gradW{2,2}.v{1} = @(u,gradU,x,y) -1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
% MANUFACTURED SOLUTION
% bottom edge
loc = @(x,y) y+1;
n = [0 -1];
tau = [1 0];
alpha = [0 1];
beta = [1 0];
eta{1} = @(x,y) alpha(1)*P(x,y) + beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2));
BC.field{1}.bndry{1} = assembleBoundary(loc,alpha,beta,eta);
% right edge
loc = @(x,y) x-1;
n = [1 0];
tau = [0 1];
alpha = [0 1];
beta = [1 0];
eta{1} = @(x,y) alpha(1)*P(x,y) + beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2));
BC.field{1}.bndry{2} = assembleBoundary(loc,alpha,beta,eta);
% top edge
loc = @(x,y) y-1;
n = [0 1];
tau = [-1 0];
alpha = [0 1];
beta = [1 0];
eta{1} = @(x,y) alpha(1)*P(x,y) + beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2));
BC.field{1}.bndry{3} = assembleBoundary(loc,alpha,beta,eta);
% left edge
loc = @(x,y) x+1;
n = [-1 0];
tau = [0 -1];
alpha = [0 1];
beta = [1 0];
eta{1} = @(x,y) alpha(1)*P(x,y) + beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2));
BC.field{1}.bndry{4} = assembleBoundary(loc,alpha,beta,eta);

%% SOLVER
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 3;

%% NULLSPACE
% pick fields that have a constant nullSpace
% nullSpace.field{2} = 'true';

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