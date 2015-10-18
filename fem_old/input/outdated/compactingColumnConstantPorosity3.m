% COMPACTINGCOLUMNCONSTANTPOROSITY.m
%
% 1D compacting column, in non-dimensionalized form, on [-L 0].  The
% multidimensional mechanical equations in this case reduce to the 1D
% equations:
%
% 1. vr = -phi*q_f'
% 2. v_r' + 1/(1-phi)*(q_f-q) = 0
% 3. q' - 4/3*(1-phi)*v_m'' = -(1-phi)
% 4. v_m' - phi/(1-phi)*(q_f-q) = 0
%
% with boundary conditions:
%
%             v_r = 0
%             v_m = 0
%               ___
%                |
%                |
%                |
%                |
%               _|_
%
%             v_r = 0
%             v_m = 0
%
% To replicate this problem in 2D we solve the full system on the 
% rectangular domain [0 1]x[-L 0]:
%
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
L = 3;
nx = 100;
% x = 0:.1:L; x = [x (0:.1:L)/10+(L-L/10)/2];  x = [x (0:.1:L)/100+(L-L/100)/2]; x = [x (0:.1:L)/1000+(L-L/1000)/2]; x = unique(sort(x));
% x = 0:.1:L; x = [x (0:.1:L)/10+(L-L/10)/2];  x = unique(sort(x));
y = [0 1];
% uniform rectangular mesh on [0 1]x[0 L]
% [nodes, edges, cells] = QuadMesh(x,y);
[nodes, edges, cells] = uniformQuadMesh(nx,1,[0 L],[0 1]);

%% VARIABLES---------------------------------------------------------------
field(1).name = 'relative velocity';
field(1).numComp = 2;
field(2).name = 'fluid pressure potential';
field(2).numComp = 1;
field(3).name = 'matrix velocity';
field(3).numComp = 2;
field(4).name = 'matrix pressure potential';
field(4).numComp = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 5;

%% FINITE ELEMENT BASIS----------------------------------------------------
field(1).basisType = @RT0Velocity;
field(2).basisType = @piecewiseConst;
field(3).basisType = @AWVelocity;
field(4).basisType = @piecewiseConst;

%% EXACT SOLUTION
phi = .02;
n = 2.2;

R = ((1-phi)*((4/3)*phi+1)*phi^(n-1))^(-1/2);
% apply boundary condition
B = exp(-R*L);
C = (1-B)/(1-B^2);

% v_r = -(1-phi)*phi^(n-1)*(C1*exp(z/(R*phi^((n-1)/2))) + C2*exp(-z/(R*phi^((n-1)/2))) + 1)

field(1).uExact{1} = @(x,y) -(1-phi)*phi^(n/2)*(1 - C*exp(-R*x) - C*exp(R*(x-L)));
field(1).uExact{2} = @(x,y) 0;
field(2).uExact{1} = @(x,y) (1-phi)*(x - (C/R)*(1 - exp(-R*L)) + (C/R)*(exp(-R*x) - exp(R*(x-L))));
field(3).uExact{1} = @(x,y) (1-phi)*phi^n*(1 - C*exp(-R*x) - C*exp(R*(x-L)));
field(3).uExact{2} = @(x,y) 0;
field(4).uExact{1} = @(x,y) (1-phi)*(x - (C/R)*(1 - exp(-R*L)) + (4*phi/(4*phi+3))*(C/R)*(exp(-R*x) - exp(R*(x-L))));


field = symbolicDerivatives(field,1);

%% VARIATIONAL FORM--------------------------------------------------------
% ( u_r , v_r ) - ( phi^(n/2)*q_f , divV_r )
field(1).varForm.v{1} = @(u, gradU, x, y) u.field{1}(1);
field(1).varForm.v{2} = @(u, gradU, x, y) u.field{1}(2);
field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -phi^(n/2)*u.field{2}(1);
field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -phi^(n/2)*u.field{2}(1);

% -( phi^(n/2)*divU_r, w_f ) - ( (phi/(1-phi))*(q_f-q) , w_f )
field(2).varForm.v{1} = @(u, gradU, x, y)  -phi^(n/2)*(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - (phi/(1-phi))*(u.field{2}(1)-u.field{4}(1));

% -( q , divV_m ) + ( 2*(1-phi)*DU_m , DV_m ) - ( (2/3)*(1-phi)*divU_m , divV_m ) - ([(1-phi) 0] , v_m)
field(3).varForm.gradV{1,1} = @(u, gradU, x, y) 2*(1-phi)*gradU.field{3}(1,1) - (2/3)*(1-phi)*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
field(3).varForm.gradV{1,2} = @(u, gradU, x, y) (1-phi)*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
field(3).varForm.gradV{2,1} = @(u, gradU, x, y) (1-phi)*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
field(3).varForm.gradV{2,2} = @(u, gradU, x, y) 2*(1-phi)*gradU.field{3}(2,2) - (2/3)*(1-phi)*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
field(3).varForm.v{1} = @(u, gradU, x, y) -(1-phi);

% -( divU_m, w ) + ( (phi/(1-phi))*(q_f-q) , w )
field(4).varForm.v{1} = @(u, gradU, x, y) -gradU.field{3}(1,1) - gradU.field{3}(2,2) + (phi/(1-phi))*(u.field{2}(1)-u.field{4}(1));

% n = 2;
% % A = .04;
% B = appCtx.B;
% C = L/2;
% D = appCtx.D;
% % phi = @(x,y) A*exp(-(x-C)^2/D^2) + B;
% % dxPhi = @(x,y) -(2*A/(D^2))*(x-C)*exp(-(x-C)^2/D);
% dyPhi = @(x,y) 0;
% 
% % phi = @(x,y) (B + (x-C)^2*(x>=0))*(sign(x-C)*(x-C)^2 <= .05 - B) + D*(sign(x-C)*(x-C)^2 > .05 - B);
%         phi = @(x,y) (B + .5*(x-1/2)^2*(x>=1/2))*(sign(x-1/2)*.5*(x-1/2)^2 <= .05 - B) + D*(sign(x-1/2)*.5*(x-1/2)^2 > .05 - B);
%         dxPhi = @(x,y) 1*(x-1/2)*(x>=1/2)*(sign(x-1/2)*.5*(x-1/2)^2 <= .05 - B);
% % dxPhi = @(x,y) 2*x*(x>=0)*(sign(x-C)*(x-L/2)^2 <= .05 - B);
% 
% 
% % phi = @(x,y) .01*(x-L/2)^2*(x>=3/2) + .00001;
% % dxPhi = @(x,y) .01*2*(x-L/2)*(x>=3/2);
% 
% % ( u_r , v_r ) - ( phi^(n/2)*q_f , divV_r )
% field(1).varForm.v{1} = @(u, gradU, x, y) u.field{1}(1) - (n/2)*phi(x,y)^((n-2)/2)*dxPhi(x,y);
% field(1).varForm.v{2} = @(u, gradU, x, y) u.field{1}(2) - (n/2)*phi(x,y)^((n-2)/2)*dyPhi(x,y);
% field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -phi(x,y)^(n/2)*u.field{2}(1);
% field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -phi(x,y)^(n/2)*u.field{2}(1);
% 
% % -( phi^(n/2)*divU_r, w_f ) - ( (phi/(1-phi))*(q_f-q) , w_f )
% field(2).varForm.v{1} = @(u, gradU, x, y)  -phi(x,y)^(n/2)*(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - ...
%     (n/2)*phi(x,y)^((n-2)/2)*(dxPhi(x,y)*u.field{1}(1) + dyPhi(x,y)*u.field{1}(2)) - ...
%     (phi(x,y)/(1-phi(x,y)))*(u.field{2}(1)-u.field{4}(1));
% 
% % -( q , divV_m ) + ( 2*(1-phi)*DU_m , DV_m ) - ( (2/3)*(1-phi)*divU_m , divV_m ) - ([(1-phi) 0] , v_m)
% field(3).varForm.gradV{1,1} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(1,1) - (2/3)*(1-phi(x,y))*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
% field(3).varForm.gradV{1,2} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
% field(3).varForm.gradV{2,1} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
% field(3).varForm.gradV{2,2} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(2,2) - (2/3)*(1-phi(x,y))*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
% field(3).varForm.v{1} = @(u, gradU, x, y) -(1-phi(x,y));
% 
% % -( divU_m, w ) + ( (phi/(1-phi))*(q_f-q) , w )
% field(4).varForm.v{1} = @(u, gradU, x, y) -(gradU.field{3}(1,1) + gradU.field{3}(2,2)) + (phi(x,y)/(1-phi(x,y)))*(u.field{2}(1)-u.field{4}(1));


%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
% EXACT SOLUTION
% U_r
% bottom edge
field(1).bndry(1).loc = @(x,y) y-1;
field(1).bndry(1).alpha = [0 1];
field(1).bndry(1).beta = [1 0];
field(1).bndry(1).eta = {@(x,y) 0 , @(x,y) 0 };
% right edge
field(1).bndry(2).loc = @(x,y) x-L;
field(1).bndry(2).alpha = [0 1];
field(1).bndry(2).beta = [1 0];
field(1).bndry(2).eta = {@(x,y) 0 , @(x,y) 0};
% top edge
field(1).bndry(3).loc = @(x,y) y;
field(1).bndry(3).alpha = [0 1];
field(1).bndry(3).beta = [1 0];
field(1).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
% left edge
field(1).bndry(4).loc = @(x,y) x;
field(1).bndry(4).alpha = [0 1];
field(1).bndry(4).beta = [1 0];
field(1).bndry(4).eta = {@(x,y) 0 , @(x,y) 0};

% U_m
% bottom edge
field(3).bndry(1).loc = @(x,y) y-1;
field(3).bndry(1).alpha = [0 1];
field(3).bndry(1).beta = [1 0];
field(3).bndry(1).eta = {@(x,y) 0 , @(x,y) 0};
% right edge
field(3).bndry(2).loc = @(x,y) x-L;
field(3).bndry(2).alpha = [0 0];
field(3).bndry(2).beta = [1 1];
field(3).bndry(2).eta = {@(x,y) 0 , @(x,y) 0};
% top edge
field(3).bndry(3).loc = @(x,y) y;
field(3).bndry(3).alpha = [0 1];
field(3).bndry(3).beta = [1 0];
field(3).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
% left edge
field(3).bndry(4).loc = @(x,y) x;
field(3).bndry(4).alpha = [0 0];
field(3).bndry(4).beta = [1 1];
field(3).bndry(4).eta = {@(x,y) 0 , @(x,y) 0};

% % bottom edge
% loc = @(x,y) y-L;
% alpha = [0 0];
% beta = [1 1];
% eta{1} = @(x,y) 0;
% eta{2} = @(x,y) 0;
% BC.field{3}.bndry{1} = assembleBoundary(loc,alpha,beta,eta);
% % right edge
% loc = @(x,y) x-1;
% alpha = [0 1];
% beta = [1 0];
% eta{1} = @(x,y) 0;
% eta{2} = @(x,y) 0;
% BC.field{3}.bndry{2} = assembleBoundary(loc,alpha,beta,eta);
% % top edge
% loc = @(x,y) y;
% alpha = [0 0];
% beta = [1 1];
% eta{1} = @(x,y) 0;
% eta{2} = @(x,y) 0;
% BC.field{3}.bndry{3} = assembleBoundary(loc,alpha,beta,eta);
% % left edge
% loc = @(x,y) x;
% alpha = [0 1];
% beta = [1 0];
% eta{1} = @(x,y) 0;
% eta{2} = @(x,y) 0;
% BC.field{3}.bndry{4} = assembleBoundary(loc,alpha,beta,eta);

%% NEWTONSOLVER------------------------------------------------------------
relTol = 1e-10;
absTol = 1e-10;
maxIter = 2;

%% NULLSPACE---------------------------------------------------------------
field(2).nullSpace = 'true';
% field(4).nullSpace = 'true';

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