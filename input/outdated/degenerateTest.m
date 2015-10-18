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
% uniform rectangular mesh on [-1 1]x[-1 1]
% nx = 10;
% ny = 10;

[nodes, edges, cells] = uniformQuadMesh(nx,ny,[-1 1],[-1 1]);

%% VARIABLES---------------------------------------------------------------
field(1).name = 'velocity';
field(1).numComp = 2;
field(2).name = 'scaled pressure';
field(2).numComp = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 5;

%% FINITE ELEMENT BASIS----------------------------------------------------
field(1).basisType = @RT0Velocity;
field(2).basisType = @piecewiseConst;

%% EXACT SOLUTION----------------------------------------------------------

visc = 1;
n = 2;
SETPHI = 3; % alpha = 2;
SETPRES = 1; % beta = 1;
SETSOURCE = 0;
switch SETPHI
    case 1
        phi = @(x,y) 1;
        dxPhi = @(x,y) 0;
        dyPhi = @(x,y) 0;
    case 2
        phi = @(x,y) (1-x^2)*(1-y^2);
        dxPhi = @(x,y) -2*x*(1-y^2);
        dyPhi = @(x,y) -(1-x^2)*2*y;
    case 3
        phi = @(x,y) x^2*y^2;
        dxPhi = @(x,y) 2*x*y^2;
        dyPhi = @(x,y) 2*x^2*y;
    case 4
        phi = @(x,y) sqrt(x)*y^2*(x>=0)*(y>=0);
        dxPhi = @(x,y) (1/2*sqrt(x))*y^2*(x>=0)*(y>=0);
        dyPhi = @(x,y) 2*sqrt(x)*y*(x>=0)*(y>=0);
    case 5
        if ~exist('alpha', 'var')
            error('myApp:argChk','undefined parameter alpha for phi test case');
        end
        phi = @(x,y) x^alpha*(x>=0);
        dxPhi = @(x,y) alpha*x^(alpha-1)*(x>=0);
        dyPhi = @(x,y) 0;
    case 6
        phi = @(x,y) sqrt(x+3/4)*(y+3/4)^2*(x>=-3/4)*(y>=3/4);
        dxPhi = @(x,y) (1/2*sqrt(x+3/4))*(y+3/4)^2*(x>=-3/4)*(y>=3/4);
        dyPhi = @(x,y) sqrt(x+3/4)*2*(y+3/4)*(x>=-3/4)*(y>=3/4);
    case 7
        phi = @(x,y) ((x+3/4)*(y+3/4)^2)^alpha*(x>=-3/4)*(y>=-3/4);
        dxPhi = @(x,y) alpha*((x+3/4)*(y+3/4)^2)^(alpha-1)*(y+3/4)^2*(x>=-3/4)*(y>=-3/4);
        dyPhi = @(x,y) alpha*((x+3/4)*(y+3/4)^2)^(alpha-1)*(x+3/4)*2*(y+3/4)*(x>=-3/4)*(y>=-3/4);
    otherwise
        error('myApp:argChk','improper phi test case number');
end

switch SETPRES
    case 1
        p = @(x,y) x^2*cos(y);
        dxP = @(x,y) 2*x*cos(y);
        dxxP = @(x,y) 2*cos(y);
        dyP = @(x,y) -x^2*sin(y);
        dyyP = @(x,y) -x^2*cos(y);
    case 2
        p = @(x,y) (1-x^2)^2*(1-y^2)^2;
        dxP = @(x,y) -4*x*(1-x^2)*(1-y^2)^2;
        dxxP = @(x,y) -4*(1-3*x^2)*(1-y^2)^2;
        dyP = @(x,y) -(1-x^2)^2*4*y*(1-y^2);
        dyyP = @(x,y) -(1-x^2)^2*4*(1-3*y^2);
    case 3
        p = @(x,y) cos(6*x^2*y);
        dxP = @(x,y) -12*x*y*sin(6*x^2*y);
        dxxP = @(x,y) -12*y*sing(6*x^2*y) - 144*x^2*y^2*cos(6*x^2*y);
        dyP = @(x,y) -6*x^2*sin(6*x^2*y);
        dyyP = @(x,y) -36*x^4*cos(6*x^2*y);
    case 4
        if ~exist('beta', 'var')
            error('myApp:argChk','undefined parameter beta for pressure test case');
        end
        r1 = (-3+sqrt(13))/2; r2 = (-3-sqrt(13))/2;
        p = @(x,y) (beta*x^r1 - r1*x^beta)/(r1*(beta-r1)*(beta-r2))*(x>=0);
        dxP = @(x,y) beta*(x^(r1-1) - x^(beta-1))/((beta-r1)*(beta-r2));
        dxxP = @(x,y) beta*( (r1-1)*x^(r1-2) - (beta-1)*x^(beta-2))/((beta-r1)*(beta-r2));
        dyP = @(x,y) 0;
        dyyP = @(x,y) 0;
    case 5
        p = @(x,y) x^2*cos(y) + y^2*cos(x);
        dxP = @(x,y) 2*x*cos(y) - y^2*sin(x);
        dxxP = @(x,y) 2*cos(y) - y^2*cos(x);
        dyP = @(x,y) -x^2*sin(y) + 2*y*cos(x);
        dyyP = @(x,y) -x^2*cos(y) + 2*cos(x);
    case 6
        p = @(x,y) cos(6*x*y^2);
        dxP = @(x,y) -6*y^2*sin(6*y^2*x);
        dxxP = @(x,y) -36*y^4*cos(6*y^2*x);
        dyP = @(x,y) -12*x*y*sin(6*y^2*x);
        dyyP = @(x,y) -12*x*sin(6*y^2*x) - 144*x^2*y^2*cos(6*y^2*x);
    otherwise
        error('myApp:argChk','improper pressure test number');
end

if SETPRES~=0
    vx = @(x,y) -phi(x,y)^(n/2)*dxP(x,y);
    vy = @(x,y) -phi(x,y)^(n/2)*dyP(x,y);
    
    field(1).uExact{1} = @(x,y) vx(x,y);
    field(1).uExact{2} = @(x,y) vy(x,y);
    field(2).uExact{1} = @(x,y) sqrt(phi(x,y))*p(x,y);
end

switch SETSOURCE
    case 0
        source = @(x,y) - visc*phi(x,y)*(dxxP(x,y) + dyyP(x,y)) - 2*visc*(dxPhi(x,y)*dxP(x,y) + dyPhi(x,y)*dyP(x,y)) + p(x,y);
    case 1
        source = @(x,y) -y;
    case 2
        source = @(x,y) x*y;
    case 3
        if ~exist('beta', 'var')
            error('myApp:argChk','undefined parameter beta for source test case');
        end
        source = @(x,y) abs(x)^beta;
    case 4
        if ~exist('beta', 'var')
            error('myApp:argChk','undefined parameter beta for source test case');
        end
        source = @(x,y) x^beta*(x>=0);
    otherwise
        error('myApp:argChk','improper source test number');
end

%% VARIATIONAL FORM--------------------------------------------------------
% (u , v) - (q , phi^(-1/2)*div(phi^(n/2)*v) =
% (u , v) - ((n/2)*phi^((n-3)/2)*gradPhi*q , v) - (phi^((n-1)/2)*q , div v)
field(1).varForm.v{1} = @(u, gradU, x, y) u.field{1}(1) -(n/2)*phi(x,y)^((n-3)/2)*dxPhi(x,y)*u.field{2}(1);
field(1).varForm.v{2} = @(u, gradU, x, y) u.field{1}(2) -(n/2)*phi(x,y)^((n-3)/2)*dyPhi(x,y)*u.field{2}(1);
field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -phi(x,y)^((n-1)/2)*u.field{2}(1);
field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -phi(x,y)^((n-1)/2)*u.field{2}(1);

% (phi^(-1/2)*div(phi^(n/2)*u , w) + (q , w) - (phi^(1/2)*source , w) =
% ((n/2)*phi^((n-3)/2)*gradPhi*u , w) - (phi^((n-1)/2)*div u , w) + (q , w) - (phi^(1/2)*source , w)
field(2).varForm.v{1} = @(u, gradU, x, y) (n/2)*phi(x,y)^((n-3)/2)*(dxPhi(x,y)*u.field{1}(1) + dyPhi(x,y)*u.field{1}(2)) + ...
    phi(x,y)^((n-1)/2)*(gradU.field{1}(1,1) + gradU.field{1}(2,2)) + u.field{2}(1) - phi(x,y)^(1/2)*source(x,y)/visc;

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0.1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
field(1).bndry(1).loc = @(x,y) (x-1)*(y-1)*(x+1)*(y+1);
field(1).bndry(1).alpha = [1 1];
field(1).bndry(1).beta = [0 0];
field(1).bndry(1).eta = {@(x,y) phi(x,y)^(n/2)*p(x,y), @(x,y) 0};

%% SOLVER------------------------------------------------------------------
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 1;

%% NULLSPACE---------------------------------------------------------------
% pick fields that have a constant nullSpace
% field(2).nullSpace = 'true';