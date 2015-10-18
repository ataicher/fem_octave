% StokesSymmetricCornerFlow.m
%
% Solve the PDE on [0,1]x[-1,0]
%   -div gradU + gradP = f
%   divU = g
%
%                               u*n = 0
%                             u*tau = -v_0 
%                        ____________________
%                       |                    |
%                       |                    |   
%                       |                    |   
%        u*n = 0        |                    | (gradU-pI)n = 0
%     tau'*gradU*n = 0  |                    | 
%                       |                    | 
%                       |                    | 
%                       |____________________|
%
%                          (gradU-pI)n = 0
%
% V = {v in H^1 : v|(y=0) = 0, v*tau|(x=0) = 0}
% Q = {q in L^2 : (q , 1) = const}
% for all v in V 
%    (gradU , gradV) - (f , v) - <n'*(gradU-pI)*n , v*n> - <tau'*gradU*n , v*tau> = 0
% for all q in Q
%    -(divU , q) + (g ,q ) = 0


% 1. (u_r , v_r) - (phi*q_f , div(phi*v_r) + <phi*q_f , v_r*n> = 0
% 
%
% 2. (div(phi*U_r , w_f) + ( phi/(1-phi)*(q_f-q) , w_f) = 0
%
% 3. -(q , divV_m) + (2*(1-phi)DU_m , DV_m) - (2/3*(1-phi)*divU_m , divV_m)
%   + <n'*(-2*(1-phi)*DU_m + (2/3*(1-phi)*divU_m+q)*I)*n , v_m*n> 
%   + <tau'*(-2*(1-phi)*DU_m)*n , v_m*tau> = 0
% 
% 4. (divU_m  , w) - (phi/(1-phi)*(q_f-q) , w) = 0
%
% for constant porosity phi, equations reduce to:
%
% 1. (u_r , v_r) - (phi*q_f , divV_r) + <phi*q_f , v_r*n> = 0
% 
% 2. (phi*divU_r , w_f) + ( phi/(1-phi)*(q_f-q) , w_f) = 0
%
% 3. -(q , divV_m) + (2*(1-phi)*DU_m , DV_m) - (2/3*(1-phi)*divU_m , divV_m)
%   + <n'*(-2*(1-phi)*DU_m + (2/3*(1-phi)*divU_m+q)*I)*n , v_m*n> 
%   + <tau'*(-2*(1-phi)*DU_m)*n , v_m*tau> = 0
% 
% 4. (divU_m  , w) - (phi/(1-phi)*(q_f-q) , w) = 0
%
% with boundary condtions:
%
%                      u_r = 0
%                      u_m = 0
%                   _____________
%                  |             |
%                  |             |
%                  |             |
%      u_r*n = 0   |             |   u_r*n = 0
%      u_m*n = 0   |             |   u_m*n = 0
%   tau'*u_m*n = 0 |             | tau'*u_m*n = 0
%                  |             |
%                  |             |
%                  |_____________|
%
%                      u_r = 0
%                      u_m = 0

%% MESH--------------------------------------------------------------------
Lx = 4;
Ly = 4;
% uniform rectangular mesh on [0 1]x[0 1]
% [nodes, edges, cells] = uniformQuadMesh(8,8,[0 1],[-1 0]);
% rectangular mesh on [0 1]x[-1 0] refined at top-left corner
nx = 10;
ny = 10;
x = Lx*(0:1/nx:1);
y = Ly*(0:1/ny:1);
% nx = Lx*[.01:.01:.1 .2:.1:.999 1];
% ny = Ly*[.01:.01:.1 .2:.1:.999 1];
% x = Lx*[0 .025 .05 .1 .15 .25 .35 .45 .55 .65 .75 .85 .95 .999 1];
% y = Ly*[0 .025 .05 .1 .15 .25 .35 .45 .55 .65 .75 .85 .95 .999 1];
% nx = Lx*[0 .025 .05 .1 .999 1];
% ny = Ly*[0 .025 .05 .1 .999 1];
% ny = [-1 -.75 -.55 -.35 -.025 0];
[nodes, edges, cells] = QuadMesh(x,y);

%% VARIABLES---------------------------------------------------------------
field(1).name = 'relative velocity';
field(1).numComp = 2;
field(2).name = 'fluid pressure potential';
field(2).numComp = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 5;

%% FINITE ELEMENT BASIS----------------------------------------------------

field(1).basisType = @RT0Velocity;
field(2).basisType = @piecewiseConst;

%% EXACT SOLUTION----------------------------------------------------------

phi = .024;
w_0 = phi^2*4.9e-6;
U_0 = 3.18e-9;
w_0/U_0

A = 0;
B = 2/pi;

r = @(x,y) sqrt(x^2+y^2);
r2 = @(x,y) x^2+y^2;
theta = @(x,y) atan(y/x);

% ux_r = -(w_0*(r^2 + 2*B*cos(2*theta)))/(U_0*r^2)
% uy_r = -(2*B*w_0*sin(2*theta))/(U_0*r^2)
% q_f = cos(theta)*(r - (2*B)/r)
field(1).uExact{1} = @(x,y) -(w_0/U_0)*(1 + 2*B*cos(2*theta(x,y))/r2(x,y));
field(1).uExact{2} = @(x,y) -(w_0/U_0)*2*B*sin(2*theta(x,y))/r2(x,y);
field(2).uExact{1} = @(x,y) cos(theta(x,y))*(r(x,y) - (2*B)/r(x,y));

field = symbolicDerivatives(field,1);

%% VARIATIONAL FORM--------------------------------------------------------
% ((w_0/w_0)*u_r , v_r) - (q_f , divV_r) = 0
field(1).varForm.v{1} = @(u, gradU, x, y) (U_0/w_0)*u.field{1}(1);
field(1).varForm.v{2} = @(u, gradU, x, y) (U_0/w_0)*u.field{1}(2);
field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -u.field{2}(1);
field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -u.field{2}(1);

% -(divU_m , w) = 0
field(2).varForm.v{1} = @(u, gradU, x, y) -gradU.field{1}(1,1) - gradU.field{1}(2,2);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0.1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
% v_r
% top and right edge
field(1).bndry(1).loc = @(x,y) (x-Lx)*(y-Ly);
field(1).bndry(1).alpha = [1 1];
field(1).bndry(1).beta = [0 0];
field(1).bndry(1).eta = {@(x,y) x , @(x,y) 0};
% bottom edge
field(1).bndry(2).loc = @(x,y) y;
field(1).bndry(2).alpha = [0 1];
field(1).bndry(2).beta = [1 0];
field(1).bndry(2).eta = {@(x,y) 0 , @(x,y) 0};
% left edge
field(1).bndry(3).loc = @(x,y) x;
field(1).bndry(3).alpha = [1 1];
field(1).bndry(3).beta = [0 0];
field(1).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};

%% SOLVER------------------------------------------------------------------
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 1;

%% NULLSPACE---------------------------------------------------------------
% pick fields that have a constant nullSpace
% field(2).nullSpace = 'true';

%% DISPLAY AND DEBUG OPTIONS-----------------------------------------------
% displayMesh = 'true';
% displayInitialGuess = 'true';
% outputResidualError = 'true';
% displayResidualError = 'true';
% displayCurrentNewtonIterate = 'true';
displayComputedSolution = 'true';
displayExactSolution = 'true';
displayStreamLines = 'true';
% displayJacobianSpyPlot = 'true';