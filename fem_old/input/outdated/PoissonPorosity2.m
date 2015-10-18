% STOKES
% Solve the PDE [-1,1]^2
%   -div gradU + gradP = f
%   -divU = -g
% for all v in V
%    (gradU , gradV) - (f , v) - <gradU.*n , v> = 0
% for all q in Q
%    (divU , q) -(g ,q ) = 0

%% MESH--------------------------------------------------------------------
% uniform rectangular mesh on [-1 1]^2
[nodes, edges, cells] = uniformQuadMesh(nx,ny,[-1 1],[-1 1]);

%% VARIABLES---------------------------------------------------------------
name.field{1} = 'velocity';
numComp.field{1} = 2;
name.field{2} = 'intermediate velocity';
numComp.field{2} = 2;
name.field{3} = 'pressure';
numComp.field{3} = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 5;

%% FINITE ELEMENT BASIS----------------------------------------------------
% basisType.field{1} = @P1Quad2;
% basisType.field{2} = @piecewiseConst;
% basisType.field{1} = @AWVelocity;
% basisType.field{2} = @AWVelocity;
% basisType.field{3} = @piecewiseConst;
basisType.field{1} = @RT0Velocity;
basisType.field{2} = @RT0Velocity;
basisType.field{3} = @piecewiseConst;

%% MANUFACTURED SOLUTION
syms x y;

% porosity
phi = @(x,y) ((y+.5)^2*(y>-.5));
gradPhix = @(x,y) 0;
gradPhiy = @(x,y) (2*(y+.5)*(y>-.5));
gradPhixx = @(x,y) 0;
gradPhiyy = @(x,y) (2*(y>-.5));

% phi = @(x,y) 1;
% gradPhix = @(x,y) 0;
% gradPhiy = @(x,y) 0;
% gradPhixx = @(x,y) 0;
% gradPhiyy = @(x,y) 0;

% manufactured pressure
% P = @(x,y) y^3*(1-y^2);
P = @(x,y) (1-x^2)*(1-y^2) + 1;
dxP = matlabFunction(diff(P,x),'vars',[x,y]);
dyP = matlabFunction(diff(P,y),'vars',[x,y]);
dxxP = matlabFunction(diff(dxP,x),'vars',[x,y]);
dxyP = matlabFunction(diff(dxP,y),'vars',[x,y]);
dyyP = matlabFunction(diff(dyP,y),'vars',[x,y]);

% obtain velocity from pressure
U2x = @(x,y) -dxP(x,y);
U2y = @(x,y) -dyP(x,y);
dxU2x = @(x,y) -dxxP(x,y);
dxU2y = @(x,y) -dxyP(x,y);
dyU2x = @(x,y) -dxyP(x,y);
dyU2y = @(x,y) -dyyP(x,y);

U1x = @(x,y) phi(x,y)^2*U2x(x,y);
U1y = @(x,y) phi(x,y)^2*U2y(x,y);
dxU1x = @(x,y) 2*phi(x,y)*gradPhix(x,y)*U2x(x,y) + phi(x,y)^2*dxU2x(x,y);
dxU1y = @(x,y) 2*phi(x,y)*gradPhix(x,y)*U2y(x,y) + phi(x,y)^2*dxU2y(x,y);
dyU1x = @(x,y) 2*phi(x,y)*gradPhiy(x,y)*U2x(x,y) + phi(x,y)^2*dyU2x(x,y);
dyU1y = @(x,y) 2*phi(x,y)*gradPhiy(x,y)*U2y(x,y) + phi(x,y)^2*dyU2y(x,y);

% set exact solution
uE.field{1} = {U1x , U1y};
gradUE.field{1} = {dxU1x , dyU1x ; dxU1y, dyU1y}; 
uE.field{2} = {U2x , U2y};
gradUE.field{2} = {dxU2x , dyU2x ; dxU2y, dyU2y}; 
uE.field{3} = {P};
gradUE.field{3} = {dxP , dyP}; 

%% VARIATIONAL FORM--------------------------------------------------------
f = @(x,y) P(x,y) + (dxU1x(x,y) + dyU1y(x,y));

% (u2 , v1) - (p , divV1)
varForm.field{1}.v{1} = @(u, gradU, x, y) u.field{2}(1);
varForm.field{1}.v{2} = @(u, gradU, x, y) u.field{2}(2);
varForm.field{1}.gradV{1,1} = @(u, gradU, x, y) -u.field{3};
varForm.field{1}.gradV{2,2} = @(u, gradU, x, y) -u.field{3};

% (u1 , v2) - (phi^2*u2 , v2)
varForm.field{2}.v{1} = @(u, gradU, x, y) u.field{1}(1) - phi(x,y)^2*u.field{2}(1);
varForm.field{2}.v{2} = @(u, gradU, x, y) u.field{1}(2) - phi(x,y)^2*u.field{2}(2);

% - (p , w) - (divU1 , w) + (f,w)
varForm.field{3}.v{1} = @(u, gradU, x, y) - u.field{3} - (gradU.field{1}(1,1) + gradU.field{1}(2,2)) + f(x,y);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 1;

% (w2 , v1) - (w3, divV1)
JacForm.field{2,1}.w{1}.v{1} = @(u,gradU,x,y) 1;
JacForm.field{2,1}.w{2}.v{2} = @(u,gradU,x,y) 1;

JacForm.field{3,1}.w{1}.v{1} = @(u,gradU,x,y) -1;
JacForm.field{3,1}.w{1}.v{2} = @(u,gradU,x,y) -1;

% (w1 , v2) - (phi^2w2, v2)
JacForm.field{1,2}.w{1}.v{1} = @(u,gradU,x,y) 1;
JacForm.field{1,2}.w{2}.v{2} = @(u,gradU,x,y) 1;

JacForm.field{2,2}.w{1}.v{1} = @(u,gradU,x,y) -phi(x,y)^2;
JacForm.field{2,2}.w{1}.v{2} = @(u,gradU,x,y) -phi(x,y)^2;

% -(w3 , v3) - (divW1 , w3)
JacForm.field{3,3}.w{1}.v{1} = @(u,gradU,x,y) -1;

JacForm.field{1,3}.gradW{1,1}.v{1} = @(u,gradU,x,y) -1;
JacForm.field{1,3}.gradW{2,2}.v{1} = @(u,gradU,x,y) -1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
% MANUFACTURED SOLUTION
loc = @(x,y) (1-x^2)*(1-y^2);
alpha = [1 1];
beta = [0 0];
eta{1} = @(x,y) alpha(1)*P(x,y);
eta{2} = @(x,y) 0;
BC.field{1}.bndry{1} = assembleBoundary(loc,alpha,beta,eta);

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