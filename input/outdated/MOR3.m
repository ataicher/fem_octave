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

%% PROBELM PARAMETERS
Theta = 0;
k_0 = 2e-10;
g = 9.8;
drho = 500;
mu_f = 1;
w_0 = k_0*drho*g/mu_f;
% w_0 = 1e-8
U_0 = 1e-9;

% mesh
%   3. uniform mesh with nx x ny elements with corner [0 x_0] x [0 y_0]
%       removed
LX = 2;
LY = 2;
nx = 16;
ny = 16;
x_0 = LX/nx;
y_0 = LY/ny;

% set porosity
%   1. constant porosity phi = phi_0
phi_0 = .04;
w_0 = w_0*phi_0^(2*(1+Theta));

%boundary conditions
%   1. Exact solution with corner removed

% naive implemenatation
%   1. divide by phi^(2+2*Theta)

%% VARIABLES---------------------------------------------------------------
field(1).name = 'Darcy velocity';
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

%% MESH--------------------------------------------------------------------
loc = @(x,y) (x >= x_0) + (y >= y_0);
Lx = [0 LX];
Ly = [0 LY];

%% POROSITY
phi = @(x,y) phi_0;
dxPhi = @(x,y) 0;
dyPhi = @(x,y) 0;

%% EXACT SOLUTION----------------------------------------------------------
% CONSTANT POROSITY
C2 = -(w_0/U_0)*phi_0^(2+2*Theta);

% u_rx = -(w_0*phi^(2+2*Theta)/U_0)*(1 + 2*B*cos(2*theta)/(r^2))
% u_ry = -(2*B*w_0*sin(2*theta))/(U_0*r^2)
field(1).uExact{1} = @(x,y) C2*(1 - 4/(pi*(x^2+y^2)^2) * (x^2-y^2) );
field(1).gradUExact{1,1} = @(x,y) 8*C2/(pi*(x^2+y^2)^3) * (-x^3 + 3*x*y^2);
field(1).gradUExact{1,2} = @(x,y) 8*C2/(pi*(x^2+y^2)^3) * (3*x^2*y - y^3);
field(1).uExact{2} = @(x,y) -4*C2/(pi*(x^2+y^2)^2) * 2*x*y ;
field(1).gradUExact{2,1} = @(x,y) 8*C2/(pi*(x^2+y^2)^3) * (3*x^2*y - y^3);
field(1).gradUExact{2,2} = @(x,y) 8*C2/(pi*(x^2+y^2)^3) * (x^3 - 3*x*y^2);

% q_f = cos(theta)*(r - (2*B)/r)
field(2).uExact{1} = @(x,y) (1-phi_0) * x*(1 - 4/(pi*(x^2+y^2)) );
field(2).gradUExact{1,1} = @(x,y) (1-phi_0) * (1 - 4/(pi*(x^2+y^2)^2) * (x^2-y^2));
field(2).gradUExact{1,2} = @(x,y) (1-phi_0) * -4/(pi*(x^2+y^2)) * 2*x*y;

% u_mx = -B*cos(theta)^2
% u_my = B*(theta - sin(theta)*cos(theta))
field(3).uExact{1} = @(x,y) -2/(pi*(x^2+y^2)) * x^2;
field(3).gradUExact{1,1} = @(x,y) -4/((x^2+y^2)^2) * x*y^2;
field(3).gradUExact{1,2} = @(x,y) 4/((x^2+y^2)^2) * -x^2*y;
field(3).uExact{2} = @(x,y) -2/(pi*(x^2+y^2)) * x*y + 2\pi * atan(y/x);
field(3).gradUExact{2,1} = @(x,y) -4/((x^2+y^2)^2) * y^3;
field(3).gradUExact{2,2} = @(x,y) -4/((x^2+y^2)^2) * x*y^2;

% q_m = cos(theta)*(r - (2*B)/r)
field(4).uExact{1} = @(x,y) (1-phi_0) * x * (1 - 4/(pi*(x^2+y^2)) );
field(4).gradUExact{1,1} = @(x,y) (1-phi_0) * (1 - 4/(pi*(x^2+y^2)) * (x^2-y^2) );
field(4).gradUExact{1,2} = @(x,y) (1-phi_0) * -4/(pi*(x^2+y^2)) * 2*x*y;

Urx = field(1).uExact{1}; Ury = field(1).uExact{2};
dxUrx = field(1).gradUExact{1,1}; dyUry = field(1).gradUExact{2,2};
Qf = field(2).uExact{1};
Umx = field(3).uExact{1}; Umy = field(3).uExact{2};
dxUmx = field(3).gradUExact{1,1}; dyUmx = field(3).gradUExact{1,2}; dxUmy = field(3).gradUExact{2,1}; dyUmy = field(3).gradUExact{2,2};
Qm = field(4).uExact{1};

%% VARIATIONAL FORM--------------------------------------------------------
% ( phi^(-2-2*Theta)*u_r , v_r ) - ( q_f , div v_r )
field(1).varForm.v{1} = @(u, gradU, x, y) (U_0/w_0)*phi(x,y)^(-2-2*Theta)*u.field{1}(1);
field(1).varForm.v{2} = @(u, gradU, x, y) (U_0/w_0)*phi(x,y)^(-2-2*Theta)*u.field{1}(2);
field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -u.field{2}(1);
field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -u.field{2}(1);

% -( div u_r, w_f ) - ( (phi/(1-phi))*(q_f-q) , w_f )
field(2).varForm.v{1} = @(u, gradU, x, y)  -(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - ...
    (phi(x,y)/(1-phi(x,y)))*(u.field{2}(1)-u.field{4}(1));

% -( q , divV_m ) + ( 2*(1-phi)*DU_m , DV_m ) - ( (2/3)*(1-phi)*divU_m , divV_m ) - ([(1-phi) ; 0] , v_m)
field(3).varForm.gradV{1,1} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(1,1) - (2/3)*(1-phi(x,y))*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
field(3).varForm.gradV{1,2} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
field(3).varForm.gradV{2,1} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
field(3).varForm.gradV{2,2} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(2,2) - (2/3)*(1-phi(x,y))*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
field(3).varForm.v{1} = @(u, gradU, x, y) -(1-phi(x,y));

% -( divU_m, w ) + ( (phi/(1-phi))*(q_f-q) , w )
field(4).varForm.v{1} = @(u, gradU, x, y) -(gradU.field{3}(1,1) + gradU.field{3}(2,2)) + (phi(x,y)/(1-phi(x,y)))*(u.field{2}(1)-u.field{4}(1));
% end

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
%EXACT SOLUTION WITH CORNER REMOVED
% v_r
% right edge
%   Dirichlet - exact solution
field(1).bndry(1).loc = @(x,y) (x-LX);
field(1).bndry(1).alpha = [0 1];
field(1).bndry(1).beta = [1 0];
field(1).bndry(1).eta = {@(x,y)  Urx(x,y), @(x,y) 0};
% top edge
%   Dirichlet - exact solution
field(1).bndry(2).loc = @(x,y) (y-LY);
field(1).bndry(2).alpha = [0 1];
field(1).bndry(2).beta = [1 0];
field(1).bndry(2).eta = {@(x,y) Ury(x,y) , @(x,y) 0};
% bottom edge
%   Dirichlet - symmetry
field(1).bndry(3).loc = @(x,y) y;
field(1).bndry(3).alpha = [0 1];
field(1).bndry(3).beta = [1 0];
field(1).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
% left edge
%   Dirichlet - exact solution
field(1).bndry(4).loc = @(x,y) x;
field(1).bndry(4).alpha = [0 1];
field(1).bndry(4).beta = [1 0];
field(1).bndry(4).eta = {@(x,y) -Urx(x,y) , @(x,y) 0};
% corner edges
%   Dirichlet - exact solution
field(1).bndry(5).loc = @(x,y) (x-x_0);
field(1).bndry(5).alpha = [0 1];
field(1).bndry(5).beta = [1 0];
field(1).bndry(5).eta = {@(x,y) -Urx(x,y), @(x,y) 0};
field(1).bndry(6).loc = @(x,y) (y-y_0);
field(1).bndry(6).alpha = [0 1];
field(1).bndry(6).beta = [1 0];
field(1).bndry(6).eta = {@(x,y) -Ury(x,y), @(x,y) 0};

% v_m
% right edge
%   Dirichlet - exact solution
field(3).bndry(1).loc = @(x,y) (x-LX);
field(3).bndry(1).alpha = [0 0];
field(3).bndry(1).beta = [1 1];
field(3).bndry(1).eta = {@(x,y) Umx(x,y) , @(x,y) Umy(x,y)};
% top edge
%   Dirichlet - exact solution
field(3).bndry(2).loc = @(x,y) (y-LY);
field(3).bndry(2).alpha = [0 0];
field(3).bndry(2).beta = [1 1];
field(3).bndry(2).eta = {@(x,y) Umy(x,y) , @(x,y) -Umx(x,y)};
% bottom edge
%   symmetry
field(3).bndry(3).loc = @(x,y) y;
field(3).bndry(3).alpha = [0 1];
field(3).bndry(3).beta = [1 0];
field(3).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
% left edge
%   Dirichlet
field(3).bndry(4).loc = @(x,y) x;
field(3).bndry(4).alpha = [0 0];
field(3).bndry(4).beta = [1 1];
field(3).bndry(4).eta = {@(x,y) 0, @(x,y) -1};
% corner edges
%   Dirichlet - exact solution
field(3).bndry(5).loc = @(x,y) (x-x_0);
field(3).bndry(5).alpha = [0 0];
field(3).bndry(5).beta = [1 1];
field(3).bndry(5).eta = {@(x,y) -Umx(x,y), @(x,y) -Umy(x,y)};
field(3).bndry(6).loc = @(x,y) (y-y_0);
field(3).bndry(6).alpha = [0 0];
field(3).bndry(6).beta = [1 1];
field(3).bndry(6).eta = {@(x,y) -Umy(x,y), @(x,y) Umx(x,y)};

%% SOLVER------------------------------------------------------------------
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 1;

%% NULLSPACE---------------------------------------------------------------
% pick fields that have a constant nullSpace
field(2).nullSpace = 'true';

%% STREAMFUNCTION
existStreamFun = 'true';
