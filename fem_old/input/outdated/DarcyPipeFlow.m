% DARCY
%   Solve the PDE on [0 L]x[-1 1]
%       u + grad p = 0
%            div u = 0
%
%                 u*n = 0
%          ____________________
%         |                    |
%         |                    |  
%   p = 1 |                    |  p = 0
%         |                    | 
%         |____________________|
%
%                 u*n = 0
%
%   for all v in V
%       (u , v) - (p , divV) + <p , v.*n>  = 0
%   for all q in Q
%       -(divU , q) = 0
%

%% MESH--------------------------------------------------------------------
L = 4;
% uniform quad mesh on a rectangular domain
[nodes, edges, cells] = uniformQuadMesh(2,3,[0 L],[-1 1]);

%% VARIABLES---------------------------------------------------------------
field(1).name = 'velocity';
field(1).numComp = 2;
field(2).name = 'pressure';
field(2).numComp = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 2;

%% FINITE ELEMENT BASIS----------------------------------------------------
field(1).basisType = @RT0Velocity;
field(2).basisType = @piecewiseConst;

%% EXACT SOLUTION
field(1).uExact{1} = @(x,y) 1/L;
field(1).uExact{2} = @(x,y) 0;
field(2).uExact{1} = @(x,y) (L-x)/L;
field = symbolicDerivatives(field);

%% VARIATIONAL FORM--------------------------------------------------------
% (u , v) - (f , v) 
field(1).varForm.v{1} = @(u, gradU, x, y) u.field{1}(1);
field(1).varForm.v{2} = @(u, gradU, x, y) u.field{1}(2);
% -(p , divV)
field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -u.field{2}(1);
field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -u.field{2}(1);
% -(divU, w) + (g , w)
field(2).varForm.v{1} = @(u, gradU, x, y) -gradU.field{1}(1,1) - gradU.field{1}(2,2);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0;
% -(w,v)
field(1).JacForm.field(1).v{1}.w{1} = @(u,gradU,x,y) 1;
field(1).JacForm.field(1).v{2}.w{2} = @(u,gradU,x,y) 1;
% -(w,divV)
field(1).JacForm.field(2).gradV{1,1}.w{1} = @(u,gradU,x,y) -1;
field(1).JacForm.field(2).gradV{2,2}.w{1} = @(u,gradU,x,y) -1;
% -(divW,v)
field(2).JacForm.field(1).v{1}.gradW{1,1} = @(u,gradU,x,y) -1;
field(2).JacForm.field(1).v{1}.gradW{2,2} = @(u,gradU,x,y) -1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
% right edge
field(1).bndry(1).loc = @(x,y) x-L;
field(1).bndry(1).alpha = [1 1];
field(1).bndry(1).beta = [0 0];
field(1).bndry(1).eta = {@(x,y) 0 , @(x,y) 0};

% top edge + bottom edge
field(1).bndry(2).loc = @(x,y) y^2-1;
field(1).bndry(2).alpha = [0 1];
field(1).bndry(2).beta = [1 0];
field(1).bndry(2).eta = {@(x,y) 0 , @(x,y) 0};

% left edge
field(1).bndry(3).loc = @(x,y) x;
field(1).bndry(3).alpha = [1 1];
field(1).bndry(3).beta = [0 0];
field(1).bndry(3).eta = {@(x,y) 1 , @(x,y) 0};

%% SOLVER
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 1;

%% NULLSPACE
% pick fields that have a constant nullSpace
% field(2).nullSpace = 'true';

%% DISPLAY AND DEBUG OPTIONS
displayMesh = 'true';
displayInitialGuess = 'true';
outputResidualError = 'true';
displayCurrentNewtonIterate = 'true';
displayComputedSolution = 'true';
displayExactSolution = 'true';
displayStreamLines = 'true';
displayJacobianSpyPlot = 'true';