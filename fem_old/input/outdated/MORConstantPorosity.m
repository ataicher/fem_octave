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
n = 2;
k_0 = 1e-9;
g = 9.8;
drho = 500;
mu_f = 1;
w_0 = k_0*drho*g/mu_f
U_0 = 1e-9
FORMULATION = 1;
Theta = 0;

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

%% MESH--------------------------------------------------------------------
SETMESH = 1;
% define mesh on [0 LX] x [0 LY]
LX = 4;
LY = 4;
%   1. uniform mesh with nx x ny elements
nx = 16;
ny = 16;
%   2. mesh refined at bottom right corner
%   3. uniform mesh with nx x ny elements with corner [0 x_0] x [0 y_0]
%       removed
x_0 = 1/32;
y_0 = 1/32;

switch SETMESH
    
    % UNIFORM MESH
    case 1
        Lx = [0 LX];
        Ly = [0 LY];
        
        % MESH REFINED AT TOP CORNER
    case 2
        x = LX*[0 .012 .025 .038 .05 .075 .1 .15 .2 .25 .35 .45 .55 .65 .75 .85 .95 1];
        y = LY*[0 .012 .025 .038 .05 .075 .1 .15 .2 .25 .35 .45 .55 .65 .75 .85 .95 1];
        
        % UNIFORM MESH WITH CORNER [0 x_0] x [0 y_0] REMOVED
    case 3
        loc = @(x,y) (x >= x_0) + (y >= y_0);
        Lx = [0 LX];
        Ly = [0 LY];
end

%% POROSITY
phi = .04;

%% EXACT SOLUTION----------------------------------------------------------

B = 2/pi;
phi_max = .04;
r = @(x,y) sqrt(x^2+y^2);
r2 = @(x,y) x^2+y^2;
theta = @(x,y) atan(y/x);

% u_rx = -(w_0*(r^2 + 2*B*cos(2*theta)))/(U_0*r^2)
% u_ry = -(2*B*w_0*sin(2*theta))/(U_0*r^2)
field(1).uExact{1} = @(x,y) -(w_0/U_0)*phi_max^(n/2)*(1 + 2*B*cos(2*theta(x,y))/r2(x,y));
field(1).uExact{2} = @(x,y) -(w_0/U_0)*phi_max^(n/2)*2*B*sin(2*theta(x,y))/r2(x,y);
% q_f = cos(theta)*(r - (2*B)/r)
field(2).uExact{1} = @(x,y) cos(theta(x,y))*(r(x,y) - (2*B)/r(x,y));
% u_mx = -B*cos(theta)^2
% u_my = B*(theta - sin(theta)*cos(theta))
field(3).uExact{1} = @(x,y) -B*cos(theta(x,y))^2;
field(3).uExact{2} = @(x,y) B*(theta(x,y) - sin(2*theta(x,y))/2);
% q_m = cos(theta)*(r - (2*B)/r)
field(4).uExact{1} = @(x,y) cos(theta(x,y))*(r(x,y) - (2*B)/r(x,y));

% field = symbolicDerivatives(field,1);

Urx = field(1).uExact{1}; Ury = field(1).uExact{2};
% dxUrx = field(1).gradUExact{1,1}; dyUry = field(1).gradUExact{2,2};
Qf = field(2).uExact{1};
Umx = field(3).uExact{1}; Umy = field(3).uExact{2};
% dxUmx = field(3).gradUExact{1,1}; dyUmx = field(3).gradUExact{1,2}; dxUmy = field(3).gradUExact{2,1}; dyUmy = field(3).gradUExact{2,2};
Qm = field(4).uExact{1};

% sigma11 = @(x,y) 2*dxUmx(x,y)-2/3*(dxUmx(x,y)+dyUmy(x,y))-Qm(x,y);
% sigma12 = @(x,y) dyUmx(x,y)+dxUmy(x,y);
% sigma21 = @(x,y) dyUmx(x,y)+dxUmy(x,y);
% sigma22 = @(x,y) 2*dyUmy(x,y)-2/3*(dxUmx(x,y)+dyUmy(x,y))-Qm(x,y);
% sigma = @(x,y) [ sigma11(x.y) sigma12(x,y) ; sigma21(x,y) sigma22(x,y)];
% t=0:.01:Lx;
% % sigma = @(x,y) dxUmx(x,y) - 2/3*(dxUmx(x,y)-dyUmy(x,y));
% sigma = @(x,y) [ 2*dxUmx(x,y)-2/3*(dxUmx(x,y)-dyUmy(x,y))-Qm(x,y)  -dyUmx(x,y)-dxUmy(x,y) ; ...
%     -dyUmx(x,y)-dxUmy(x,y) 2*dyUmy(x,y)-2/3*(dxUmx(x,y)-dyUmy(x,y))-Qm(x,y)];
% sigmaVec = zeros(size(t));
% for i=1:length(t)
%     A = sigma(t(i),Ly);
%     n = [0 1]';
%     tan = [-1 0]';
% sigmaVecN(i) = n'*A*n + t(i);
% sigmaVecTan(i) = tan'*A*n;
% end
% figure
% plot(t,sigmaVecN)
% figure
% plot(t,sigmaVecTan)

%% VARIATIONAL FORM--------------------------------------------------------
% ( u_r , v_r ) - ( phi^(n/2)*q_f , divV_r )
field(1).varForm.v{1} = @(u, gradU, x, y) (U_0/w_0)*u.field{1}(1);
field(1).varForm.v{2} = @(u, gradU, x, y) (U_0/w_0)*u.field{1}(2);
field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -phi^(n/2)*u.field{2}(1);
field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -phi^(n/2)*u.field{2}(1);

% -( phi^(n/2)*divU_r, w_f ) - ( phi*(q_f-q) , w_f )
field(2).varForm.v{1} = @(u, gradU, x, y)  -phi^(n/2)*(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - phi*(u.field{2}(1)-u.field{4}(1));

% -( q , divV_m ) + ( 2*DU_m , DV_m ) - ( (2/3)*divU_m , divV_m ) - ([1 0] , v_m)
field(3).varForm.gradV{1,1} = @(u, gradU, x, y) 2*(1-phi)*gradU.field{3}(1,1) - (2/3)*(1-phi)*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
field(3).varForm.gradV{1,2} = @(u, gradU, x, y) (1-phi)*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
field(3).varForm.gradV{2,1} = @(u, gradU, x, y) (1-phi)*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
field(3).varForm.gradV{2,2} = @(u, gradU, x, y) 2*(1-phi)*gradU.field{3}(2,2) - (2/3)*(1-phi)*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
field(3).varForm.v{1} = @(u, gradU, x, y) -1;

% -( divU_m, w ) + ( (phi/(1-phi))*(q_f-q) , w )
field(4).varForm.v{1} = @(u, gradU, x, y) -gradU.field{3}(1,1) - gradU.field{3}(2,2) + phi*(u.field{2}(1)-u.field{4}(1));

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0.1;

%% SOLVER------------------------------------------------------------------
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 1;

%% NULLSPACE---------------------------------------------------------------
% pick fields that have a constant nullSpace
field(2).nullSpace = 'true';
field(4).nullSpace = 'true';

%% BOUNDARY CONDITIONS-----------------------------------------------------
SETBC = 1;
% boundary condition options
%   1. Exact solution with corner removed
%   2. Exact solution
%   3. zero stress far field condtion with corner removed
%   4. zero stress far field condition
%   5. zero stress far field condition with no fluid flow through left
%   boundary except at the bottom [0 y_f]
y_f = 4/16;

% run boundary condition script
MORBC;