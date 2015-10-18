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

%% MESH--------------------------------------------------------------------
SETMESH = 3;
% mesh on [0 1]^2
%   1. uniform mesh with nx x ny elements
nx = 16;
ny = 16;
%   2. uniform mesh with nx x ny elements with corner [0 x_0] x [0 y_0]
%       removed
x_0 = 1/32;
y_0 = 1/32;
%   3. mesh refined at bottom right corner
switch SETMESH
    
    % UNIFORM MESH [0 1]^2
    case 1
        Lx = [0 1];
        Ly = [0 1];
    
    % UNIFORM MESH ON [0 1]^2 \ [0 1/8]^2
    case 2
        Lx = [0 1];
        Ly = [0 1];
        loc = @(x,y) (x >= x_0) + (y >= y_0); 
        
    % MESH ON [0 1]^2 REFINED AT BOTTOM LEFT CORNER
    case 3
        x = [0 .025 .05 .1 .15 .25 .35 .45 .55 .65 .75 .85 .95 1];
        y = [0 .025 .05 .1 .15 .25 .35 .45 .55 .65 .75 .85 .95 1];
end

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
% field(1).basisType = @TaylorHoodVelocity;
% field(2).basisType = @P1Quad;

%% EXACT SOLUTION----------------------------------------------------------
B = 2/pi;
r = @(x,y) sqrt(x^2+y^2);
r2 = @(x,y) x^2+y^2;
theta = @(x,y) atan(y/x);

% u_x = B*cos(theta)^2
field(1).uExact{1} = @(x,y) -B*cos(theta(x,y))^2;
field(1).gradUExact{1,1} = @(x,y) -2*B*y*cos(theta(x,y))*sin(theta(x,y))/r2(x,y);
field(1).gradUExact{1,2} = @(x,y) 2*B*x*cos(theta(x,y))*sin(theta(x,y))/r2(x,y);
% u_y = B*(theta - sin(theta)*cos(theta))
field(1).uExact{2} = @(x,y) B*(theta(x,y) - sin(2*theta(x,y))/2);
field(1).gradUExact{2,1} = @(x,y) -B*y*(1-cos(2*theta(x,y)))/r2(x,y);
field(1).gradUExact{2,2} = @(x,y) B*x*(1-cos(2*theta(x,y)))/r2(x,y);
% p = -2*B*U_0*cos(theta)/r
field(2).uExact{1} = @(x,y) -2*B*cos(theta(x,y))/r(x,y);
field(2).gradUExact{1,1} = @(x,y) 0;
field(2).gradUExact{1,2} = @(x,y) 0;

Ux = field(1).uExact{1}; Uy = field(1).uExact{2};
P = field(2).uExact{1};

%% VARIATIONAL FORM--------------------------------------------------------
% (gradU , gradV) - (p , divV)
field(1).varForm.gradV{1,1} = @(u, gradU, x, y) gradU.field{1}(1,1) - u.field{2}(1);
field(1).varForm.gradV{1,2} = @(u, gradU, x, y) gradU.field{1}(1,2);
field(1).varForm.gradV{2,1} = @(u, gradU, x, y) gradU.field{1}(2,1);
field(1).varForm.gradV{2,2} = @(u, gradU, x, y) gradU.field{1}(2,2) - u.field{2}(1);
% -(divU , v)
field(2).varForm.v{1} = @(u, gradU, x, y) -gradU.field{1}(1,1) - gradU.field{1}(2,2);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0.1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
SETBC = 2;

switch SETBC
    
    % EXACT SOLUTION WITH CORNER REMOVED
    case 1
        % top edge
        field(1).bndry(1).loc = @(x,y) y-1;
        field(1).bndry(1).alpha = [0 0];
        field(1).bndry(1).beta = [1 1];
        field(1).bndry(1).eta = {@(x,y) Uy(x,y) , @(x,y) -Ux(x,y)};
        % right edge
        field(1).bndry(2).loc = @(x,y) x-1;
        field(1).bndry(2).alpha = [0 0];
        field(1).bndry(2).beta = [1 1];
        field(1).bndry(2).eta = {@(x,y) Ux(x,y) , @(x,y) Uy(x,y)};
        % bottom edge
        field(1).bndry(3).loc = @(x,y) y;
        field(1).bndry(3).alpha = [0 1];
        field(1).bndry(3).beta = [1 0];
        field(1).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        field(1).bndry(4).loc = @(x,y) x;
        field(1).bndry(4).alpha = [0 0];
        field(1).bndry(4).beta = [1 1];
        field(1).bndry(4).eta = {@(x,y) 0 , @(x,y) -1};
        % corner edges
        field(1).bndry(5).loc = @(x,y) (x-x_0);
        field(1).bndry(5).alpha = [0 0];
        field(1).bndry(5).beta = [1 1];
        field(1).bndry(5).eta = {@(x,y) -Ux(x,y), @(x,y) -Uy(x,y)};
        field(1).bndry(6).loc = @(x,y) (y-y_0);
        field(1).bndry(6).alpha = [0 0];
        field(1).bndry(6).beta = [1 1];
        field(1).bndry(6).eta = {@(x,y) -Uy(x,y), @(x,y) Ux(x,y)};
        
    % EXACT SOLUTION
    case 2
        % top edge
        field(1).bndry(1).loc = @(x,y) y-1;
        field(1).bndry(1).alpha = [0 0];
        field(1).bndry(1).beta = [1 1];
        field(1).bndry(1).eta = {@(x,y) Uy(x,y) , @(x,y) -Ux(x,y)};
        % right edge
        field(1).bndry(2).loc = @(x,y) x-1;
        field(1).bndry(2).alpha = [0 0];
        field(1).bndry(2).beta = [1 1];
        field(1).bndry(2).eta = {@(x,y) Ux(x,y) , @(x,y) Uy(x,y)};
        % bottom edge
        field(1).bndry(3).loc = @(x,y) y;
        field(1).bndry(3).alpha = [0 1];
        field(1).bndry(3).beta = [1 0];
        field(1).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        field(1).bndry(4).loc = @(x,y) x;
        field(1).bndry(4).alpha = [0 0];
        field(1).bndry(4).beta = [1 1];
        field(1).bndry(4).eta = {@(x,y) 0, @(x,y) -1};
        
    % NO STRESS FAR FIELD CONDITION
    case 3
        % top and right edge
        field(1).bndry(1).loc = @(x,y) (x-1)*(y-1);
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
        field(1).bndry(3).alpha = [0 0];
        field(1).bndry(3).beta = [1 1];
        field(1).bndry(3).eta = {@(x,y) 0, @(x,y) -1};
        
    % NO STRESS FAR FIELD CONDITION WITH CORNER REMOVED
    case 4
        % top and right edge
        field(1).bndry(1).loc = @(x,y) (x-1)*(y-1);
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
        field(1).bndry(3).alpha = [0 0];
        field(1).bndry(3).beta = [1 1];
        field(1).bndry(3).eta = {@(x,y) 0, @(x,y) -1};
        % corner edges
        field(1).bndry(4).loc = @(x,y) (x-1/8);
        field(1).bndry(4).alpha = [0 0];
        field(1).bndry(4).beta = [1 1];
        field(1).bndry(4).eta = {@(x,y) -Ux(x,y), @(x,y) -Uy(x,y)};
        field(1).bndry(5).loc = @(x,y) (y-1/8);
        field(1).bndry(5).alpha = [0 0];
        field(1).bndry(5).beta = [1 1];
        field(1).bndry(5).eta = {@(x,y) -Uy(x,y), @(x,y) Ux(x,y)};
end

%% SOLVER------------------------------------------------------------------
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 5;

%% NULLSPACE---------------------------------------------------------------
% pick fields that have a constant nullSpace
field(2).nullSpace = 'true';