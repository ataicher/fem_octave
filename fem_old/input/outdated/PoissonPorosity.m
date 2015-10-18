% STOKES
% Solve the PDE [-1,1]^2
%   u + gradP = f
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
name.field{2} = 'pressure';
numComp.field{2} = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 5;

%% FINITE ELEMENT BASIS----------------------------------------------------
basisType.field{1} = @RT0Velocity;
basisType.field{2} = @piecewiseConst;
% basisType.field{1} = @AWVelocity;
% basisType.field{2} = @piecewiseConst;
% basisType.field{1} = @TaylorHoodVelocity;
% basisType.field{2} = @P1Quad;

%% MANUFACTURED SOLUTION---------------------------------------------------
syms x y;

% porosity
phi = @(x,y) (y-.5)^2*(y<.5);
gradPhix = @(x,y) 0;
gradPhiy = @(x,y) 2*(y-.5)*(y<.5);
gradPhixx = @(x,y) 0;
gradPhiyy = @(x,y) 2*(y<.5);

phi = @(x,y) y^2;
gradPhix = @(x,y) 0;
gradPhiy = @(x,y) 2*y;
gradPhixx = @(x,y) 0;
gradPhiyy = @(x,y) 2;

phi = @(x,y) 1;
gradPhix = @(x,y) 0;
gradPhiy = @(x,y) 0;
gradPhixx = @(x,y) 0;
gradPhiyy = @(x,y) 0;

% manufactured pressure
% P = @(x,y) y^3*(1-y^2);
% P = @(x,y) 1+y;
% P = @(x,y) cos(pi*x)*cos(pi*y);
P = @(x,y) 1;
dxP = matlabFunction(diff(P,x),'vars',[x,y]);
dyP = matlabFunction(diff(P,y),'vars',[x,y]);
dxxP = matlabFunction(diff(dxP,x),'vars',[x,y]);
dxyP = matlabFunction(diff(dxP,y),'vars',[x,y]);
dyyP = matlabFunction(diff(dyP,y),'vars',[x,y]);

% obtain velocity from pressure
Ux = @(x,y) -phi(x,y)*dxP(x,y);
Uy = @(x,y) -phi(x,y)*dyP(x,y);
dxUx = @(x,y) -(gradPhix(x,y)*dxP(x,y) + phi(x,y)*dxxP(x,y));
dxUy = @(x,y) -(gradPhix(x,y)*dyP(x,y) + phi(x,y)*dxyP(x,y));
dyUx = @(x,y) -(gradPhiy(x,y)*dxP(x,y) + phi(x,y)*dxyP(x,y));
dyUy = @(x,y) -(gradPhiy(x,y)*dyP(x,y) + phi(x,y)*dyyP(x,y));

% set exact solution
uE.field{1} = {Ux , Uy};
gradUE.field{1} = {dxUx , dyUx ; dxUy, dyUy}; 
uE.field{2} = {P};
gradUE.field{2} = {dxP , dyP}; 

%% VARIATIONAL FORM--------------------------------------------------------
f = @(x,y) P(x,y) + phi(x,y)*(dxUx(x,y) + dyUy(x,y)) + (gradPhix(x,y)*Ux(x,y) + gradPhiy(x,y)*Uy(x,y));

% (u , v) - (p , div(phi*v)) = (u , v) - (p*gradPhi , v) - (p*phi, divV)
varForm.field{1}.v{1} = @(u, gradU, x, y) u.field{1}(1) - u.field{2}*gradPhix(x,y);
varForm.field{1}.v{2} = @(u, gradU, x, y) u.field{1}(2) - u.field{2}*gradPhiy(x,y);
varForm.field{1}.gradV{1,1} = @(u, gradU, x, y) -u.field{2}*phi(x,y);
varForm.field{1}.gradV{2,2} = @(u, gradU, x, y) -u.field{2}*phi(x,y);

% - (p , w) - (div(phi*u) , w) + (f,w)
varForm.field{2}.v{1} = @(u, gradU, x, y) - u.field{2} - phi(x,y)*(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - (gradPhix(x,y)*u.field{1}(1) + gradPhiy(x,y)*u.field{1}(2)) + f(x,y);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0;

% (w1 , v1) - (w2*gradPhi , v1) - (w2*phi, divV1)
JacForm.field{1,1}.w{1}.v{1} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.w{2}.v{2} = @(u,gradU,x,y) 1;

JacForm.field{2,1}.w{1}.v{1} = @(u,gradU,x,y) -gradPhix(x,y);
JacForm.field{2,1}.w{1}.v{2} = @(u,gradU,x,y) -gradPhiy(x,y);

JacForm.field{2,1}.w{1}.gradV{1,1} = @(u,gradU,x,y) -phi(x,y);
JacForm.field{2,1}.w{1}.gradV{2,2} = @(u,gradU,x,y) -phi(x,y);

% -(w2 , v2) - (phi*divW1 , w2) - (gradPhi.*w1 , w2)
JacForm.field{2,2}.w{1}.v{1} = @(u,gradU,x,y) -1;

JacForm.field{1,2}.gradW{1,1}.v{1} = @(u,gradU,x,y) -phi(x,y);
JacForm.field{1,2}.gradW{2,2}.v{1} = @(u,gradU,x,y) -phi(x,y);

JacForm.field{1,2}.w{1}.v{1} = @(u,gradU,x,y) -gradPhix(x,y);
JacForm.field{1,2}.w{2}.v{1} = @(u,gradU,x,y) -gradPhiy(x,y);

%% BOUNDARY CONDITIONS-----------------------------------------------------
% MANUFACTURED SOLUTION
% bottom edge
loc = @(x,y) y+1;
alpha = [1 1];
beta = [0 0];
eta{1} = @(x,y) phi(x,y)*P(x,y);
eta{2} = @(x,y) 0;
BC.field{1}.bndry{1} = assembleBoundary(loc,alpha,beta,eta);

% right edge
loc = @(x,y) x-1;
alpha = [1 1];
beta = [0 0];
eta{1} = @(x,y) phi(x,y)*P(x,y);
eta{2} = @(x,y) 0;
BC.field{1}.bndry{2} = assembleBoundary(loc,alpha,beta,eta);

% top edge
loc = @(x,y) y-1;
alpha = [1 1];
beta = [0 0];
eta{1} = @(x,y) phi(x,y)*P(x,y);
eta{2} = @(x,y) 0;
BC.field{1}.bndry{3} = assembleBoundary(loc,alpha,beta,eta);

% left edge
loc = @(x,y) x+1;
alpha = [1 1];
beta = [0 0];
eta{1} = @(x,y) phi(x,y)*P(x,y);
eta{2} = @(x,y) 0;
BC.field{1}.bndry{4} = assembleBoundary(loc,alpha,beta,eta);

%% SOLVER------------------------------------------------------------------
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 3;

%% NULLSPACE---------------------------------------------------------------
% pick fields that have a constant nullSpace
% nullSpace.field{2} = 'true';

%% DISPLAY AND DEBUG OPTIONS-----------------------------------------------
% displayMesh = 'true';
% displayInitialGuess = 'true';
% displayResidualError = 'true';
% displayCurrentNewtonIterate = 'true';
displayComputedSolution = 'true';
exactSolutionExist = 'true';
displayExactSolution = 'true';
% displayStreamLines = 'true';
% displayJacobianSpyPlot = 'true';