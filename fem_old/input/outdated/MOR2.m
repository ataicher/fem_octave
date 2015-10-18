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
k_0 = 1e-9;
g = 9.8;
drho = 500;
mu_f = 1;
w_0 = k_0*drho*g/mu_f;
U_0 = 1e-9;

% mesh
SETMESH = 1;
% define mesh on [0 LX] x [0 LY]
%   1. uniform mesh with nx x ny elements
%   2. mesh refined at bottom right corner
%   3. uniform mesh with nx x ny elements with corner [0 x_0] x [0 y_0]
%       removed
LX = 2;
LY = 2;
if ~exist('nx','var') || ~exist('ny','var')
    nx = 10;
    ny = 10;
end
x_0 = LX/nx;
y_0 = LY/ny;

% set porosity
if ~exist('SETPHI','var')
    SETPHI = 1;
end
%   1. constant porosity phi = phi_0
%   2. MOR-like porosity
if ~exist('phi_0','var')
    phi_0 = .04;
end
phi_min = .001;
phi_max = .04;

%boundary conditions
SETBC = 1;
% boundary condition options
%   1. Exact solution with corner removed
%   2. Exact solution
%   3. zero stress far field condtion with corner removed
%   4. zero stress far field condition
%   5. zero stress far field condition with no fluid flow through left
%       boundary except at the bottom [0 y_f]
y_f = 4/16;

% degenerate problem
if ~exist('SETDEGENERATE','var')
    SETDEGENERATE = 0;
end
% under-integration
SETUNDERINTEGRATE = 0;


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

%% DEGENERATE--------------------------------------------------------------
if SETDEGENERATE
    field(2).name = 'scaled fluid pres. pot.';
    degenerate = 'true';
    if SETUNDERINTEGRATE
        underIntegrate = 'true';
    end
    RT0VelocityField = 1;
    RT0PressureField = 2;
end

%% MESH--------------------------------------------------------------------
switch SETMESH
    
    % UNIFORM MESH
    case 1
        Lx = [0 LX];
        Ly = [0 LY];
        
        % MESH REFINED AT BOTTOM-LEFT CORNER
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
switch SETPHI
    % CONSTANT POROSITY
    case 1
        phi = @(x,y) phi_0;
        dxPhi = @(x,y) 0;
        dyPhi = @(x,y) 0;
        
        % MOR-LIKE POROSITY
    case 2
        x_1 = LX/nx;
        y_1 = x_1;
        x_2 = 4*LX/nx;
        
        phi = @(x,y) (phi_min - phi_max)*( (y>=y_1)*(x<=x_1) + (y>=x)*(x>=x_1)*(x<=x_2) ) + phi_max;
        dxPhi = @(x,y) 0;
        dyPhi = @(x,y) 0;
        
                % phi = @(x,y) (phi_min - phi_max)*( (y>=y_1)*(x<=x_1) + (y>=x)*(x>=x_1)*(x<=x_2) ) + phi_max;
        phi = @(x,y) (phi_min - phi_max)*( (x>=x_1)*(y<=y_1) + (x>=x_1+nx)*(y>y_1)*(y<=y_1+ny) + (x>=x_1+2*nx)*(y>y_1+ny)*(y<=y_1+2*ny) + (x>=x_1+3*nx)*(y>y_1+2*ny)*(y<=y_1+3*ny) + (x>=x_1+4*nx)*(y>y_1+3*ny)*(y<=y_1+4*ny));
        dxPhi = @(x,y) 0;
        dyPhi = @(x,y) 0;
end

%% EXACT SOLUTION----------------------------------------------------------
switch SETPHI
    
    % CONSTANT POROSITY
    case 1
    B = 2/pi;
    
    r = @(x,y) sqrt(x^2+y^2);
    r2 = @(x,y) x^2+y^2;
    theta = @(x,y) atan(y/x);
    
    % u_rx = -(w_0*(r^2 + 2*B*cos(2*theta)))/(U_0*r^2))
    field(1).uExact{1} = @(x,y) -(w_0/U_0)*phi_0^(1+Theta)*(1 + 2*B*cos(2*theta(x,y))/r2(x,y));
    field(1).gradUExact{1,1} = @(x,y) -(w_0/U_0)*phi_0^(1+Theta)*((4*B+2)*x^3 + (12*B-2)*x*y^2)/((x^2+y^2)^3);
    field(1).gradUExact{1,2} = @(x,y) -(w_0/U_0)*phi_0^(1+Theta)*((4*B-2)*y^3 - (12*B+2)*x^2*y)/((x^2+y^2)^3);
    % u_ry = -(2*B*w_0*sin(2*theta))/(U_0*r^2)
    field(1).uExact{2} = @(x,y) -(w_0/U_0)*phi_0^(1+Theta)*2*B*sin(2*theta(x,y))/r2(x,y);
    field(1).gradUExact{2,1} = @(x,y) -4*B*(w_0/U_0)*phi_0^(1+Theta)*(y^3-3*x^2*y)/((x^2+y^2)^3);
    field(1).gradUExact{2,2} = @(x,y) -4*B*(w_0/U_0)*phi_0^(1+Theta)*(x^3-3*x*y^2)/((x^2+y^2)^3);
    
    if SETDEGENERATE
        % q_f = sqrt(phi_0)*cos(theta)*(r - (2*B)/r)
        field(2).uExact{1} = @(x,y) sqrt(phi_0)*cos(theta(x,y))*(r(x,y) - (2*B)/r(x,y));
        field(2).gradUExact{1,1} = @(x,y) 0;
        field(2).gradUExact{1,2} = @(x,y) 0;
    else
        % q_f = cos(theta)*(r - (2*B)/r)
        field(2).uExact{1} = @(x,y) cos(theta(x,y))*(r(x,y) - (2*B)/r(x,y));
        field(2).gradUExact{1,1} = @(x,y) 0;
        field(2).gradUExact{1,2} = @(x,y) 0;
        % field(2).gradUExact{1,1} = @(x,y) 1 + 4*B*x/((x^2+y^2)^2);
        % field(2).gradUExact{1,1} = @(x,y) 4*B*y/((x^2+y^2)^2);
    end
    
    % u_mx = -B*cos(theta)^2
    field(3).uExact{1} = @(x,y) -B*cos(theta(x,y))^2;
    field(3).gradUExact{1,1} = @(x,y) -2*B*x*y^2/((x^2+y^2)^2);
    field(3).gradUExact{1,2} = @(x,y) -2*B*x^2*y/((x^2+y^2)^2);
    % u_my = B*(theta - sin(theta)*cos(theta))
    field(3).uExact{2} = @(x,y) B*(theta(x,y) - sin(2*theta(x,y))/2);
    field(3).gradUExact{2,1} = @(x,y) -2*B*(y^3 - 2*x^2*y)/((x^2+y^2)^2);
    field(3).gradUExact{2,2} = @(x,y) -2*B*x*y^2/((x^2+y^2)^2);
    
    % q_m = cos(theta)*(r - (2*B)/r)
    field(4).uExact{1} = @(x,y) cos(theta(x,y))*(r(x,y) - (2*B)/r(x,y));
    field(4).gradUExact{1,1} = @(x,y) 0;
    field(4).gradUExact{1,2} = @(x,y) 0;
    % field(4).gradUExact{1,1} = @(x,y) 1 + 4*B*x/((x^2+y^2)^2);
    % field(4).gradUExact{1,1} = @(x,y) 4*B*y/((x^2+y^2)^2);
    
%     field = symbolicDerivatives(field,1);
    
    Urx = field(1).uExact{1}; Ury = field(1).uExact{2};
    dxUrx = field(1).gradUExact{1,1}; dyUry = field(1).gradUExact{2,2};
    Qf = field(2).uExact{1};
    Umx = field(3).uExact{1}; Umy = field(3).uExact{2};
    dxUmx = field(3).gradUExact{1,1}; dyUmx = field(3).gradUExact{1,2}; dxUmy = field(3).gradUExact{2,1}; dyUmy = field(3).gradUExact{2,2};
    Qm = field(4).uExact{1};
end

%% VARIATIONAL FORM--------------------------------------------------------


if SETDEGENERATE
    
    % ( u_r , v_r ) - ( q_f , phi^(-1/2)*div(phi^(1+Theta)*V_r) )
    if ~SETUNDERINTEGRATE
        field(1).varForm.v{1} = @(u, gradU, x, y) (U_0/w_0)*u.field{1}(1);
        field(1).varForm.v{2} = @(u, gradU, x, y) (U_0/w_0)*u.field{1}(2);
    end
    
    % ( phi^(-1/2)*div(phi^(1+Theta)*u_r), w_f ) + ( (phi/(1-phi))*(q_f-q) , w_f )
    field(2).varForm.v{1} = @(u, gradU, x, y) (1/(1-phi(x,y)))*(u.field{2}(1) - sqrt(phi(x,y))*u.field{4}(1));
    
    % -( q , divV_m ) + ( 2*(1-phi)*DU_m , DV_m ) - ( (2/3)*(1-phi)*divU_m , divV_m ) - ([(1-phi) ; 0] , v_m)
    field(3).varForm.gradV{1,1} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(1,1) - (2/3)*(1-phi(x,y))*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
    field(3).varForm.gradV{1,2} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
    field(3).varForm.gradV{2,1} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
    field(3).varForm.gradV{2,2} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(2,2) - (2/3)*(1-phi(x,y))*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
    field(3).varForm.v{1} = @(u, gradU, x, y) -(1-phi(x,y));
    
    % -( divU_m, w ) + ( (phi/(1-phi))*(q_f-q) , w )
    field(4).varForm.v{1} = @(u, gradU, x, y) -(gradU.field{3}(1,1) + gradU.field{3}(2,2)) + (1/(1-phi(x,y)))*(sqrt(phi(x,y))*u.field{2}(1)-phi(x,y)*u.field{4}(1));
else
    
    % ( u_r , v_r ) - ( q_f , div(phi^(1+Theta)*v_r )
    field(1).varForm.v{1} = @(u, gradU, x, y) (U_0/w_0)*u.field{1}(1) - (1+Theta)*phi(x,y)^Theta*dxPhi(x,y)*u.field{2}(1);
    field(1).varForm.v{2} = @(u, gradU, x, y) (U_0/w_0)*u.field{1}(2) - (1+Theta)*phi(x,y)^Theta*dyPhi(x,y)*u.field{2}(1);
    field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -phi(x,y)^(1+Theta)*u.field{2}(1);
    field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -phi(x,y)^(1+Theta)*u.field{2}(1);
    
    % -( div(phi^(1+Theta)*u_r), w_f ) - ( (phi/(1-phi))*(q_f-q) , w_f )
    field(2).varForm.v{1} = @(u, gradU, x, y)  -phi(x,y)^(1+Theta)*(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - ...
        (1+Theta)*phi(x,y)^Theta*(dxPhi(x,y)*u.field{1}(1) + dyPhi(x,y)*u.field{1}(2)) - ...
        (phi(x,y)/(1-phi(x,y)))*(u.field{2}(1)-u.field{4}(1));
    
    % -( q , divV_m ) + ( 2*(1-phi)*DU_m , DV_m ) - ( (2/3)*(1-phi)*divU_m , divV_m ) - ([(1-phi) ; 0] , v_m)
    field(3).varForm.gradV{1,1} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(1,1) - (2/3)*(1-phi(x,y))*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
    field(3).varForm.gradV{1,2} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
    field(3).varForm.gradV{2,1} = @(u, gradU, x, y) (1-phi(x,y))*(gradU.field{3}(1,2) + gradU.field{3}(2,1));
    field(3).varForm.gradV{2,2} = @(u, gradU, x, y) 2*(1-phi(x,y))*gradU.field{3}(2,2) - (2/3)*(1-phi(x,y))*(gradU.field{3}(1,1) + gradU.field{3}(2,2)) - u.field{4}(1);
    field(3).varForm.v{1} = @(u, gradU, x, y) -(1-phi(x,y));
    
    % -( divU_m, w ) + ( (phi/(1-phi))*(q_f-q) , w )
    field(4).varForm.v{1} = @(u, gradU, x, y) -(gradU.field{3}(1,1) + gradU.field{3}(2,2)) + (phi(x,y)/(1-phi(x,y)))*(u.field{2}(1)-u.field{4}(1));
end
%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
% run boundary condition script
MORBC;

%% SOLVER------------------------------------------------------------------
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 1;

%% NULLSPACE---------------------------------------------------------------
% pick fields that have a constant nullSpace
switch SETBC
    case 1
field(2).nullSpace = 'true';
% field(4).nullSpace = 'true';
    case 2
field(2).nullSpace = 'true';
end
