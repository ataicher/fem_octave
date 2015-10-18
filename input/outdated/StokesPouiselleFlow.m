% StokesPouiselleFlow.m
%
% Solve the PDE [0 L]x[-1,1]
%   -div gradU + gradP = f
%   divU = g
%   
%                                u = 0
%                        ____________________
%                       |                    |
%      u*tau = 0        |                    |   u*tau = 0
%   n'*(gradU-pI)*n = 0 |                    |   n'*(gradU-pI)*n = 1
%                       |                    | 
%                       |____________________|
%
%                                u = 0
%   
% for all v in V
%    (gradU , gradV) - (p, divV) - (f , v) - <n'*(gradU-pI)*n , v*n> - <tau'*gradU*n , v*tau> = 0
% for all q in Q
%    -(divU , q) + (g ,q ) = 0
%

%% MESH--------------------------------------------------------------------
L = 4;
% uniform rectangular mesh on [0 L]x[-1 1]
[nodes, edges, cells] = uniformQuadMesh(nx,ny,[0 L],[-1 1]);

%% VARIABLES---------------------------------------------------------------
name.field{1} = 'velocity';
numComp.field{1} = 2;
name.field{2} = 'pressure';
numComp.field{2} = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 5;

%% FINITE ELEMENT BASIS----------------------------------------------------
% basisType.field{1} = @P1Quad2;
% basisType.field{2} = @piecewiseConst;
% basisType.field{1} = @AWVelocity;
% basisType.field{2} = @piecewiseConst;
basisType.field{1} = @TaylorHoodVelocity;
basisType.field{2} = @P1Quad;

%% MANUFACTURED SOLUTION
uE.field{1}{1} = @(x,y) (y^2-1)/(2*L);
uE.field{1}{2} = @(x,y) 0;
uE.field{2}{1} = @(x,y) (x+1)/L;

[gradUE, grad2UE] = symDerivatives(uE);
Ux = uE.field{1}{1}; Uy = uE.field{1}{2};
dxUx = gradUE.field{1}{1,1}; dyUx = gradUE.field{1}{1,2}; dxUy = gradUE.field{1}{2,1}; dyUy = gradUE.field{1}{2,2};
dxxUx = grad2UE.field{1}{1,1,1}; dyyUx = grad2UE.field{1}{1,2,2}; dxxUy = grad2UE.field{1}{2,1,1}; dyyUy = grad2UE.field{1}{2,2,2};
P = uE.field{2}{1};
dxP = gradUE.field{2}{1,1}; dyP = gradUE.field{2}{1,2};

%% VARIATIONAL FORM--------------------------------------------------------
% (gradU , gradV) - (p , divV)
varForm.field{1}.gradV{1,1} = @(u, gradU, x, y) gradU.field{1}(1,1) - u.field{2};
varForm.field{1}.gradV{1,2} = @(u, gradU, x, y) gradU.field{1}(1,2);
varForm.field{1}.gradV{2,1} = @(u, gradU, x, y) gradU.field{1}(2,1);
varForm.field{1}.gradV{2,2} = @(u, gradU, x, y) gradU.field{1}(2,2) - u.field{2};
% -(divU , v)
varForm.field{2}.v{1} = @(u, gradU, x, y) -gradU.field{1}(1,1) - gradU.field{1}(2,2);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0;

% (gradW1, gradV1)
JacForm.field{1,1}.gradW{1,1}.gradV{1,1} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.gradW{1,2}.gradV{1,2} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.gradW{2,1}.gradV{2,1} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.gradW{2,2}.gradV{2,2} = @(u,gradU,x,y) 1;

% -(w2,divV1)
JacForm.field{2,1}.w{1}.gradV{1,1} = @(u,gradU,x,y) -1;
JacForm.field{2,1}.w{1}.gradV{2,2} = @(u,gradU,x,y) -1;

% -(divW1,v2)
JacForm.field{1,2}.gradW{1,1}.v{1} = @(u,gradU,x,y) -1;
JacForm.field{1,2}.gradW{2,2}.v{1} = @(u,gradU,x,y) -1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
% right edge
loc = @(x,y) x-2;
alpha = [1 0];
beta = [0 1];
eta = {@(x,y) 1, @(x,y) 0};
BC.field{1}.bndry{1} = assembleBoundary(loc,alpha,beta,eta);
% top edge + top edge
loc = @(x,y) y*(y-L);
alpha = [0 0];
beta = [1 1];
eta = {@(x,y) 0 ; @(x,y) 0};
BC.field{1}.bndry{2} = assembleBoundary(loc,alpha,beta,eta);
% left edge
loc = @(x,y) x+2;
alpha = [1 0];
beta = [0 1];
eta = {@(x,y) 0 ; @(x,y) 0};
BC.field{1}.bndry{3} = assembleBoundary(loc,alpha,beta,eta);

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