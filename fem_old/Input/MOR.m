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
w_0 = 1;
U_0 = .005;

% mesh on [0 LX] x [0 LY]
%   1. uniform mesh with nx x ny elements
%   2. mesh refined at bottom right corner
%   3. uniform mesh with nx x ny elements with corner [0 x_0] x [0 y_0]
%       removed
SETMESH = 3;
LX = 2;
LY = 2;
if ~exist('nx','var') || ~exist('ny','var')
    nx = 16;
    ny = 16;
end
x_0 = LX/nx;
y_0 = LY/ny;

% set porosity
%   1. constant porosity phi = phi_0
%   2. MOR-like porosity
if ~exist('SETPHI','var')
    SETPHI = 1;
end
phi_0 = .04;
if ~exist('phi_min','var')
    phi_min = 0;
end
phi_max = .04;

% boundary conditions
%   1. Exact solution with corner removed
%   2. Exact solution
%   3. zero stress far field condtion with corner removed
%   4. zero stress far field condition
%   5. zero stress far field condition with no fluid flow through left
%       boundary except at the bottom [0 y_f]
SETBC = 1;
y_f = LY/ny;

% pick implemenatation
%   1. standard (divide by phi^(2+2*Theta))
%   2. auxilliary velocity
%   3. symmetric
%   4. scaled
if ~exist('FORMULATION','var')
    FORMULATION = 1;
end

%% MESH--------------------------------------------------------------------
switch SETMESH
    
    % uniform mesh
    case 1
        Lx = [0 LX];
        Ly = [0 LY];
        
        % mesh refined at bottom left corner
    case 2
        x = LX*[0 .012 .025 .038 .05 .075 .1 .15 .2 .25 .35 .45 .55 .65 .75 .85 .95 1];
        y = LY*[0 .012 .025 .038 .05 .075 .1 .15 .2 .25 .35 .45 .55 .65 .75 .85 .95 1];
        
        % uniform mesh with corner [0 x_0] x [0 y_0] removed
    case 3
        loc = @(x,y) (x >= x_0) + (y >= y_0);
        Lx = [0 LX];
        Ly = [0 LY];
end

%% POROSITY----------------------------------------------------------------
switch SETPHI
    % constant
    case 1
        phi = @(x,y) phi_0;
        dxPhi = @(x,y) 0;
        dyPhi = @(x,y) 0;
        
        % MOR--like
    case 2
        x_1 = 2*LX/nx;
        y_1 = LY/ny;
        
        % phi = @(x,y) (phi_min - phi_max)*( (y>=y_1)*(x<=x_1) + (y>=x)*(x>=x_1)*(x<=x_2) ) + phi_max;
        phi = @(x,y) (phi_min - phi_max)*( (x>=x_1)*(y<=y_1) + (x>=x_1+nx)*(y>y_1)*(y<=y_1+ny) + ...
            (x>=x_1+2*nx)*(y>y_1+ny)*(y<=y_1+2*ny) + (x>=x_1+3*nx)*(y>y_1+2*ny)*(y<=y_1+3*ny) + (x>=x_1+4*nx)*(y>y_1+3*ny)*(y<=y_1+4*ny));
        dxPhi = @(x,y) 0;
        dyPhi = @(x,y) 0;
end

%% EXACT SOLUTION----------------------------------------------------------
    % constant porosity
if SETPHI == 1 && SETBC <= 4   
    MORConstantPhiExactSol
end

%% VARIATIONAL FORM--------------------------------------------------------
mantleMechanicsVarForm

%% VARIABLES---------------------------------------------------------------
if FORMULATION == 3 || 4
    field(1).name = 'scaled relative velocity';
else
    field(1).name = 'Darcy velocity';
end
field(1).numComp = 2;
if FORMULATION == 4
    field(2).name = 'scaled fluid pres. pot.';
else
    field(2).name = 'fluid pressure pot.';
end
field(2).numComp = 1;
field(3).name = 'matrix velocity';
field(3).numComp = 2;
field(4).name = 'matrix pressure pot.';
field(4).numComp = 1;
if FORMULATION == 2
    field(5).name = 'auxilliary velocity';
    field(5).numComp = 2;
end

%% GUASS QUADRATURE--------------------------------------------------------
prec = 5;

%% FINITE ELEMENT BASIS----------------------------------------------------
field(1).basisType = @RT0Velocity;
field(2).basisType = @piecewiseConst;
field(3).basisType = @AWVelocity;
field(4).basisType = @piecewiseConst;
if FORMULATION == 2
    field(5).basisType = @RT0Velocity;
end

%% DEGENERATE--------------------------------------------------------------
if FORMULATION == 4
    SETDEGENERATE = 1;
    degenerate = 'true';
    % under-integration
    if SETUNDERINTEGRATE
        underIntegrate = 'true';
    end
    RT0VelocityField = 1;
    RT0PressureField = 2;
else
    SETDEGENERATE = 0;
end

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
MORBC2;

%% SOLVER------------------------------------------------------------------
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 1;

%% NULLSPACE---------------------------------------------------------------
if SETBC == 1
    field(2).nullSpace = 'true';
    field(4).nullSpace = 'true';
end

%% STREAMFUNCTION
if SETPHI == 1
    existStreamFun = 'true';
end