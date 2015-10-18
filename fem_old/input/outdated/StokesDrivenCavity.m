% STOKESDRIVENCAVITY
%
% Solve the PDE [-1,1]^2
%   -div gradU + gradP = 0
%   divU = 0
%
%                 u*tau = -1
%           ____________________
%          |                    |
%          |                    |   
%   u = 0  |                    | u = 0
%          |                    | 
%          |____________________|
%
%                  u = 0
%
%
% for all v in V
%    (gradU , gradV) - (p, divV) - (f , v) - <n'*(gradU-pI).*n , v*n> - <tau'*gradU.*n , v*tau> = 0
% for all q in Q
%    -(divU , q)  = 0

%% MESH--------------------------------------------------------------------
% uniform rectangular mesh on [-1 1]^2
[nodes, edges, cells] = uniformQuadMesh(5,10,[0 1],[-2 0]);

%% VARIABLES---------------------------------------------------------------
field(1).name = 'velocity';
field(1).numComp = 2;
field(2).name = 'pressure';
field(2).numComp = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 5;

%% FINITE ELEMENT BASIS----------------------------------------------------
field(1).basisType = @AWVelocity;
field(2).basisType = @piecewiseConst;

%% VARIATIONAL FORM--------------------------------------------------------
% (gradU , gradV) - (p , divV)
field(1).varForm.gradV{1,1} = @(u, gradU, x, y) gradU.field{1}(1,1) - u.field{2};
field(1).varForm.gradV{1,2} = @(u, gradU, x, y) gradU.field{1}(1,2);
field(1).varForm.gradV{2,1} = @(u, gradU, x, y) gradU.field{1}(2,1);
field(1).varForm.gradV{2,2} = @(u, gradU, x, y) gradU.field{1}(2,2) - u.field{2};
% -(divU , v)
field(2).varForm.v{1} = @(u, gradU, x, y) -gradU.field{1}(1,1) - gradU.field{1}(2,2);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0.1;
% (gradW, gradV)
field(1).JacForm.field(1).gradV{1,1}.gradW{1,1} = @(u,gradU,x,y) 1;
field(1).JacForm.field(1).gradV{1,2}.gradW{1,2} = @(u,gradU,x,y) 1;
field(1).JacForm.field(1).gradV{2,1}.gradW{2,1} = @(u,gradU,x,y) 1;
field(1).JacForm.field(1).gradV{2,2}.gradW{2,2} = @(u,gradU,x,y) 1;
% -(w,divV)
field(1).JacForm.field(2).v{1}.gradW{1,1} = @(u,gradU,x,y) -1;
field(1).JacForm.field(2).v{1}.gradW{2,2} = @(u,gradU,x,y) -1;
% -(divW,v)
field(2).JacForm.field(1).gradV{1,1}.w{1} = @(u,gradU,x,y) -1;
field(2).JacForm.field(1).gradV{2,2}.w{1} = @(u,gradU,x,y) -1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
% bottom edge
field(1).bndry(1).loc = @(x,y) y+2;
field(1).bndry(1).alpha = [0 0];
field(1).bndry(1).beta = [1 1];
field(1).bndry(1).eta = {@(x,y) 0, @(x,y) 0};

% right edge
field(1).bndry(2).loc = @(x,y) x-1;
field(1).bndry(2).alpha = [0 0];
field(1).bndry(2).beta = [1 1];
field(1).bndry(2).eta = {@(x,y) 0, @(x,y) 0};

% top edge
field(1).bndry(3).loc = @(x,y) y;
field(1).bndry(3).alpha = [0 0];
field(1).bndry(3).beta = [1 1];
field(1).bndry(3).eta = {@(x,y) 0, @(x,y) 1};

% left edge
field(1).bndry(4).loc = @(x,y) x;
field(1).bndry(4).alpha = [0 0];
field(1).bndry(4).beta = [1 1];
field(1).bndry(4).eta = {@(x,y) 0, @(x,y) 0};

%% SOLVER
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 1;

%% NULLSPACE
% pick fields that have a constant nullSpace
field(2).nullSpace = 'true';

%% DISPLAY AND DEBUG OPTIONS
displayMesh = 'true';
displayInitialGuess = 'true';
outputResidualError = 'true';
displayResidualError = 'true';
displayCurrentNewtonIterate = 'true';
displayComputedSolution = 'true';
displayExactSolution = 'true';
displayStreamLines = 'true';
displayJacobianSpyPlot = 'true';