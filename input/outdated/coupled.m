% Coupled.m
%
% Solve the PDE
%   -div(gradU1) = f1
%   -div(gradU2) + (u1 - u2) = f2
% with Dirichlet boundary conditions
%
% for all v1 in V1
%    (gradU1 , gradV1) - (f1 , v1) - <gradU1.*n , v1> = 0
% for all v2 in  V2
%   (gradU2 , gradV2) + (u1-u2 ,v2) - (f2 , v2) - <gradU2.*n , v2>= 0 

%% MESH--------------------------------------------------------------------
% uniform quad mesh on a rectangular domain
% [nodes, edges, cells] = uniformQuadMesh(10,10,[-1 1],[-1 1]);
[nodes, edges, cells] = QuadMesh([-1 -.9 -.8 -.4  1],[-1 .4 .8 .9 1]);

%% VARIABLES---------------------------------------------------------------
name.field{1} = 'u1';
numComp.field{1} = 2;
name.field{1} = 'u2';
numComp.field{2} = 2;

%% GUASS QUADRATURE--------------------------------------------------------
% polynomial order of accuracy
prec = 2;

%% FINITE ELEMENT BASIS----------------------------------------------------
basisType.field{1} = @P1Quad2;
basisType.field{2} = @P1Quad2;

%% MANUFACTURED SOLUTION
uE.field{1}{1} = @(x,y) x*y;
uE.field{1}{2} = @(x,y) x*y;
uE.field{2}{1} = @(x,y) x;
uE.field{2}{2} = @(x,y) x*y;

[gradUE, grad2UE] = symDerivatives(uE);
U1x = uE.field{1}{1}; U1y = uE.field{1}{2}; 
dxU1x = gradUE.field{1}{1,1}; dyU1x = gradUE.field{1}{1,2}; dxU1y = gradUE.field{1}{2,1}; dyU1y = gradUE.field{1}{2,2};
dxxU1x = grad2UE.field{1}{1,1,1}; dyyU1x = grad2UE.field{1}{1,2,2}; dxxU1y = grad2UE.field{1}{2,1,1}; dyyU1y = grad2UE.field{1}{2,2,2};
U2x = uE.field{2}{1}; U2y = uE.field{2}{2}; 
dxU2x = gradUE.field{2}{1,1}; dyU2x = gradUE.field{2}{1,2}; dxU2y = gradUE.field{2}{2,1}; dyU2y = gradUE.field{2}{2,2};
dxxU2x = grad2UE.field{2}{1,1,1}; dyyU2x = grad2UE.field{2}{1,2,2}; dxxU2y = grad2UE.field{2}{2,1,1}; dyyU2y = grad2UE.field{2}{2,2,2};

%% VARIATIONAL FORM--------------------------------------------------------
% manufactured solution forcing functions
f1x = @(x,y) -(dxxU1x(x,y) + dyyU1x(x,y));
f1y = @(x,y) -(dxxU1y(x,y) + dyyU1y(x,y));
f2x = @(x,y) -(dxxU2x(x,y) + dyyU2x(x,y)) + (U1x(x,y) - U2x(x,y));
f2y = @(x,y) -(dxxU2y(x,y) + dyyU2y(x,y)) + (U1y(x,y) - U2y(x,y));
% -(f1 , v1)
varForm.field{1}.v{1} = @(u,gradU,x,y) -f1x(x,y);
varForm.field{1}.v{2} = @(u,gradU,x,y) -f1y(x,y);
% ((1 + mu*|u2|^2)*gradU1 , gradV1)
varForm.field{1}.gradV{1,1} = @(u,gradU,x,y) gradU.field{1}(1,1);
varForm.field{1}.gradV{1,2} = @(u,gradU,x,y) gradU.field{1}(1,2);
varForm.field{1}.gradV{2,1} = @(u,gradU,x,y) gradU.field{1}(2,1);
varForm.field{1}.gradV{2,2} = @(u,gradU,x,y) gradU.field{1}(2,2);
% (u1 - u2 - f2 , v2)
varForm.field{2}.v{1} = @(u,gradU,x,y) u.field{1}(1) - u.field{2}(1) - f2x(x,y);
varForm.field{2}.v{2} = @(u,gradU,x,y) u.field{1}(2) - u.field{2}(2) - f2y(x,y);
% (gradU2 , gradV2)
varForm.field{2}.gradV{1,1} = @(u,gradU,x,y) gradU.field{2}(1,1);
varForm.field{2}.gradV{1,2} = @(u,gradU,x,y) gradU.field{2}(1,2);
varForm.field{2}.gradV{2,1} = @(u,gradU,x,y) gradU.field{2}(2,1);
varForm.field{2}.gradV{2,2} = @(u,gradU,x,y) gradU.field{2}(2,2);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0;
% (gradW1 , gradV1)
JacForm.field{1,1}.gradW{1,1}.gradV{1,1} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.gradW{1,2}.gradV{1,2} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.gradW{2,1}.gradV{2,1} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.gradW{2,2}.gradV{2,2} = @(u,gradU,x,y) 1;
% (w1 , v2)
JacForm.field{1,2}.w{1}.v{1} = @(u,gradU,x,y) 1;
JacForm.field{1,2}.w{2}.v{2} = @(u,gradU,x,y) 1;
% -(w2 , v2)
JacForm.field{2,2}.w{1}.v{1} = @(u,gradU,x,y) -1;
JacForm.field{2,2}.w{2}.v{2} = @(u,gradU,x,y) -1;
% (gradW2 , gradV2)
JacForm.field{2,2}.gradW{1,1}.gradV{1,1} = @(u,gradU,x,y) 1;
JacForm.field{2,2}.gradW{1,2}.gradV{1,2} = @(u,gradU,x,y) 1;
JacForm.field{2,2}.gradW{2,1}.gradV{2,1} = @(u,gradU,x,y) 1;
JacForm.field{2,2}.gradW{2,2}.gradV{2,2} = @(u,gradU,x,y) 1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
alpha = [0 0];
beta = [1 1];

% bottom edge
loc = @(x,y) y+1;
n = [0 -1];
tau = [1 0];
%
eta{1} = @(x,y) beta(1)*(U1x(x,y)*n(1) + U1y(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(U1x(x,y)*tau(1) + U1y(x,y)*tau(2));
BC.field{1}.bndry{1} = assembleBoundary(loc,alpha,beta,eta);
%
eta{1} = @(x,y) beta(1)*(U2x(x,y)*n(1) + U2y(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(U2x(x,y)*tau(1) + U2y(x,y)*tau(2));
BC.field{2}.bndry{1} = assembleBoundary(loc,alpha,beta,eta);

% right edge
loc = @(x,y) x-1;
n = [1 0];
tau = [0 1];
%
eta{1} = @(x,y) beta(1)*(U1x(x,y)*n(1) + U1y(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(U1x(x,y)*tau(1) + U1y(x,y)*tau(2));
BC.field{1}.bndry{2} = assembleBoundary(loc,alpha,beta,eta);
%
eta{1} = @(x,y) beta(1)*(U2x(x,y)*n(1) + U2y(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(U2x(x,y)*tau(1) + U2y(x,y)*tau(2));
BC.field{2}.bndry{2} = assembleBoundary(loc,alpha,beta,eta);

% top edge
loc = @(x,y) y-1;
n = [0 1];
tau = [-1 0];
eta{1} = @(x,y) beta(1)*(U1x(x,y)*n(1) + U1y(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(U1x(x,y)*tau(1) + U1y(x,y)*tau(2));
BC.field{1}.bndry{3} = assembleBoundary(loc,alpha,beta,eta);
%
eta{1} = @(x,y) beta(1)*(U2x(x,y)*n(1) + U2y(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(U2x(x,y)*tau(1) + U2y(x,y)*tau(2));
BC.field{2}.bndry{3} = assembleBoundary(loc,alpha,beta,eta);

% left edge
loc = @(x,y) x+1;
n = [-1 0];
tau = [0 -1];
%
eta{1} = @(x,y) beta(1)*(U1x(x,y)*n(1) + U1y(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(U1x(x,y)*tau(1) + U1y(x,y)*tau(2));
BC.field{1}.bndry{4} = assembleBoundary(loc,alpha,beta,eta);
%
eta{1} = @(x,y) beta(1)*(U2x(x,y)*n(1) + U2y(x,y)*n(2));
eta{2} = @(x,y) beta(2)*(U2x(x,y)*tau(1) + U2y(x,y)*tau(2));
BC.field{2}.bndry{4} = assembleBoundary(loc,alpha,beta,eta);

%% SOLVER------------------------------------------------------------------
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 3;

%% NULLSPACE---------------------------------------------------------------
% pick fields that have a constant nullSpace

%% DISPLAY AND DEBUG OPTIONS-----------------------------------------------
displayMesh = 'true';
displayInitialGuess = 'true';
% displayResidualError = 'true';
% displayCurrentNewtonIterate = 'true';
displayComputedSolution = 'true';
exactSolutionExist = 'true';
displayExactSolution = 'true';
% displayStreamLines = 'true';
% displayJacobianSpyPlot = 'true';