% COMPACTINGCOLUMNCONSTANTPOROSITY.m
%
% 1D compacting column, in non-dimensionalized form, on [-L 0].  The
% multidimensional mechanical equations in this case reduce to the 1D
% equations:
%
% 1. vr = -phi^(1+Theta)*q_f'
% 2. (phi^(1+Theta)v_r)' + phi/(1-phi)*(q_f-q) = 0
% 3. (q - 4/3*(1-phi)*v_m')' = 1-phi
% 4. v_m' - phi/(1-phi)*(q_f-q) = 0
%
% with boundary conditions:
% 1.
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
%
% 2.
%       phi^(1+Theta)*q_f = 0
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
% rectangular domain [-L L]x[0 1]:
%
% 1. (u_r , v_r) - (phi*q_f , div(phi^(1*v_r)) + <phi*q_f , v_r*n> = 0
%
% 2. (div(phi*U_r) , w_f) + ( phi/(1-phi)*(q_f-q) , w_f) = 0
%
% 3. -(q , divV_m) + (2*(1-phi)DU_m , DV_m) - (2/3*(1-phi)*divU_m , divV_m)
%   + <n'*(-2*(1-phi)*DU_m + (2/3*(1-phi)*divU_m+q)*I)*n , v_m*n>
%   + <tau'*(-2*(1-phi)*DU_m)*n , v_m*tau> = 0
%
% 4. (divU_m  , w) - (phi/(1-phi)*(q_f-q) , w) = 0
%
% with boundary condtions:
% 1.
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
%
%
% 2.
%                phi^(1+Theta)*q_f = 0
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
Theta = .25;

% mesh
SETMESH = 1;
% rectangular mesh on [-L L]x[0 1]
% 1. uniform mesh with nx x 1 elements
% 2. mesh refined around x = L/2
L = .2;
if ~exist('nx','var')
    nx = 50;
end

% porosity
if ~exist('SETPHI','var')
    SETPHI = 3;
end
%   1. constant: phi = phi_0
%   2. Gaussian: phi = A*exp(-x^2/B^2) + phi_0
%   3. Quadratic: (only for Theta = 0)
%             { phi_m               x <= 0
%       phi = { phi_m + C*x^2   0 < x <= ((phi_p-phi_m)/C)^(1/2)
%             { phi_p               x > ((phi_p-phi_m)/C)^(1/2)
%   4. Square root:
%             { phi_m                   x <= 0
%       phi = { phi_m + D*x^(1/2)   0 < x <= ((phi_p-phi_m)/D)^2
%             { phi_p                   x > ((phi_p-phi_m)/D)^2
%   5. Discontinuous jump:
%       phi = {phi_m x <  0
%             {phi_p x >= 0
phi_0 = 0.04;
if ~exist('phi_m','var')
    phi_m = 0;
end
phi_p = .04;
A = .04;
B = .2;
C = 1;
D = .08;

% boundary condition
SETBC = 1;
%   1. no flow at both top and bottom boundaries
%   2. zero fluid pressure at top boundary and no flow at bottom

% pick naive implemenatation
%   1. divide by phi^(2+2*Theta)
%   2. symmetric with phi
%   3. antisymmetric
%   4. auxiliary variable
if ~exist('NAIVE','var')
    NAIVE = 1;
end

% degenerate problem
if ~exist('SETDEGENERATE','var')
    SETDEGENERATE = 1;
end
% under-integration
SETUNDERINTEGRATE = 0;

%% MESH--------------------------------------------------------------------
switch SETMESH
    
    % UNIFORM MESH
    case 1
        ny = 1;
        Lx = [-L L]; Ly = [0 1];
        
        % REFINED AT CENTER
    case 2
        x = [-L:L/50:-L/10 -L/10+L/1000:L/1000:L/10-L/1000 L/10:L/50:L];
        y = [0 1];
end

%% VARIABLES---------------------------------------------------------------
if SETDEGENERATE || NAIVE == 2
    field(1).name = 'scaled relative velocity';
else
    field(1).name = 'Darcy velocity';
    field(1).name = 'scaled relative velocity';
end
field(1).numComp = 2;
if SETDEGENERATE
    field(2).name = 'scaled fluid pres. pot.';
else
    field(2).name = 'fluid pressure pot.';
end
field(2).numComp = 1;
field(3).name = 'matrix velocity';
field(3).numComp = 2;
field(4).name = 'matrix pressure pot.';
field(4).numComp = 1;
if NAIVE == 4 && ~SETDEGENERATE
    field(5).name = 'relative velocity';
    field(5).numComp = 2;
end

%% GUASS QUADRATURE--------------------------------------------------------
prec = 5;

%% FINITE ELEMENT BASIS----------------------------------------------------
field(1).basisType = @RT0Velocity;
field(2).basisType = @piecewiseConst;
field(3).basisType = @AWVelocity;
field(4).basisType = @piecewiseConst;
if NAIVE == 4 && ~SETDEGENERATE
    field(5).basisType = @RT0Velocity;
end

%% DEGENERATE--------------------------------------------------------------
if SETDEGENERATE
    degenerate = 'true';
    if SETUNDERINTEGRATE
        underIntegrate = 'true';
    end
    RT0VelocityField = 1;
    RT0PressureField = 2;
end

%% POROSITY
switch SETPHI
    % CONSTANT
    case 1
        phi = @(x,y) phi_0;
        dxPhi = @(x,y) 0;
        
        % GAUSSIAN
    case 2
        phi = @(x,y) A*exp(-x^2/(B^2)) + phi_0;
        dxPhi = @(x,y) -(2*A*x/(B^2))*exp(-x^2/(B^2));
        
        % QUADRATIC (only for Theta = 0)
    case 3
        phi = @(x,y) phi_m*(x<0) + (phi_m + C*x^2)*(x>=0)*(x<=sqrt((phi_p-phi_m)/C))  + phi_p*(x > sqrt((phi_p-phi_m)/C));
        dxPhi = @(x,y) 2*C*x*(x>0)*(x<=sqrt((phi_p-phi_m)/C));
        
        % SQRT
    case 4
        phi = @(x,y) phi_m*(x<0) + (phi_m + D*sqrt(x))*(x>=0)*(x<=((phi_p-phi_m)/D)^2) + phi_p*(x>((phi_p-phi_m)/D)^2);
        dxPhi = @(x,y) D/(2*sqrt(x))*(x>0)*(x<=((phi_p-phi_m)/D)^2);
        
        % DISCONTINOUS JUMP
    case 5
        phi = @(x,y) phi_m*(x<0) + phi_p*(x>=0);
        dxPhi = @(x,y) 0;
end
dyPhi = @(x,y) 0;

% store in appCtx
appCtx.phi = phi;
appCtx.dxPhi = dxPhi;
appCtx.Theta = Theta;

%% EXACT SOLUTION
switch SETPHI
    % CONSTANT
    case 1
        constPhiCompactingColumnExact;
        
    % QUADRATIC APPROXIMATION
    case 3
        if C == 1
            quadPhiCompactingColumnExact;
        end
        
    % DISCONTINUOUS JUMP
    case 5
        discPhiCompactingColumnExact;
end

%% VARIATIONAL FORM--------------------------------------------------------
if SETDEGENERATE
    
    % ( u_r , v_r ) - ( phi^(n/2)*q_f , phi^(-1/2)*div(phi^(1+Theta)*V_r) )
    if ~SETUNDERINTEGRATE
        field(1).varForm.v{1} = @(u, gradU, x, y) u.field{1}(1);
        field(1).varForm.v{2} = @(u, gradU, x, y) u.field{1}(2);
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
    
    switch NAIVE
        case 1
            % ( phi^(-2-2*Theta)*u_r , v_r ) - ( q_f , div v_r )
            field(1).varForm.v{1} = @(u, gradU, x, y) phi(x,y)^(-2-2*Theta)*u.field{1}(1);
            field(1).varForm.v{2} = @(u, gradU, x, y) phi(x,y)^(-2-2*Theta)*u.field{1}(2);
            field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -u.field{2}(1);
            field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -u.field{2}(1);
            
            % -( div u_r, w_f ) - ( (phi/(1-phi))*(q_f-q) , w_f )
            field(2).varForm.v{1} = @(u, gradU, x, y)  -(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - ...
                (phi(x,y)/(1-phi(x,y)))*(u.field{2}(1)-u.field{4}(1));
        case 2
            % ( u_r , v_r ) - ( q_f , div(phi^(1+Theta)*v_r )
            field(1).varForm.v{1} = @(u, gradU, x, y) u.field{1}(1) - (1+Theta)*phi(x,y)^Theta*dxPhi(x,y)*u.field{2}(1);
            field(1).varForm.v{2} = @(u, gradU, x, y) u.field{1}(2) - (1+Theta)*phi(x,y)^Theta*dyPhi(x,y)*u.field{2}(1);
            field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -phi(x,y)^(1+Theta)*u.field{2}(1);
            field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -phi(x,y)^(1+Theta)*u.field{2}(1);
            
            % -( div(phi^(1+Theta)*u_r), w_f ) - ( (phi/(1-phi))*(q_f-q) , w_f )
            field(2).varForm.v{1} = @(u, gradU, x, y)  -phi(x,y)^(1+Theta)*(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - ...
                (1+Theta)*phi(x,y)^Theta*(dxPhi(x,y)*u.field{1}(1) + dyPhi(x,y)*u.field{1}(2)) - ...
                (phi(x,y)/(1-phi(x,y)))*(u.field{2}(1)-u.field{4}(1));
        case 3
            % ( u_r , v_r ) - ( q_f , div(phi^(2+2*Theta)*v_r )
            field(1).varForm.v{1} = @(u, gradU, x, y) u.field{1}(1) - (2+2*Theta)*phi(x,y)^(1+2*Theta)*dxPhi(x,y)*u.field{2}(1);
            field(1).varForm.v{2} = @(u, gradU, x, y) u.field{1}(2) - (2+2*Theta)*phi(x,y)^(1+2*Theta)*dyPhi(x,y)*u.field{2}(1);
            field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -phi(x,y)^(2+2*Theta)*u.field{2}(1);
            field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -phi(x,y)^(2+2*Theta)*u.field{2}(1);
            
            % -( div u_r, w_f ) - ( (phi/(1-phi))*(q_f-q) , w_f )
            field(2).varForm.v{1} = @(u, gradU, x, y)  -(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - ...
                (phi(x,y)/(1-phi(x,y)))*(u.field{2}(1)-u.field{4}(1));
        case 4
            % ( u_r , v_r ) - ( q_f , div v_r )
            field(1).varForm.v{1} = @(u, gradU, x, y) u.field{5}(1);
            field(1).varForm.v{2} = @(u, gradU, x, y) u.field{5}(2);
            field(1).varForm.gradV{1,1} = @(u, gradU, x, y) -u.field{2}(1);
            field(1).varForm.gradV{2,2} = @(u, gradU, x, y) -u.field{2}(1);
            
            % -( div u_r, w_f ) - ( (phi/(1-phi))*(q_f-q) , w_f )
            field(2).varForm.v{1} = @(u, gradU, x, y)  -(gradU.field{1}(1,1) + gradU.field{1}(2,2)) - ...
                (phi(x,y)/(1-phi(x,y)))*(u.field{2}(1)-u.field{4}(1));
            
            %(phi^(2*(1+Theta))*su_r, sv_r) - (u_r,sv_r)
            field(5).varForm.v{1} = @(u, gradU, x, y) phi(x,y)^(2*(1+Theta))*u.field{5}(1) - u.field{1}(1);
            field(5).varForm.v{2} = @(u, gradU, x, y) phi(x,y)^(2*(1+Theta))*u.field{5}(2) - u.field{1}(2);
    end
    
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
switch SETBC
    
    % NO FLOW AT TOP AND BOTTOM
    case 1
        % u_r
        % bottom edge
        field(1).bndry(1).loc = @(x,y) y-1;
        field(1).bndry(1).alpha = [0 1];
        field(1).bndry(1).beta = [1 0];
        field(1).bndry(1).eta = {@(x,y) 0 , @(x,y) 0};
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
        field(1).bndry(4).loc = @(x,y) x+L;
        field(1).bndry(4).alpha = [0 1];
        field(1).bndry(4).beta = [1 0];
        field(1).bndry(4).eta = {@(x,y) 0 , @(x,y) 0};
        
        % u_m
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
        field(3).bndry(4).loc = @(x,y) x+L;
        field(3).bndry(4).alpha = [0 0];
        field(3).bndry(4).beta = [1 1];
        field(3).bndry(4).eta = {@(x,y) 0 , @(x,y) 0};
        
        % ZERO FLUID PRESSURE AT BOTTOM AND NO FLOW AT BOTTOM
    case 2
        % u_r
        % bottom edge
        field(1).bndry(1).loc = @(x,y) y-1;
        field(1).bndry(1).alpha = [0 1];
        field(1).bndry(1).beta = [1 0];
        field(1).bndry(1).eta = {@(x,y) 0 , @(x,y) 0};
        % right edge
        field(1).bndry(2).loc = @(x,y) x-L;
        field(1).bndry(2).alpha = [1 1];
        field(1).bndry(2).beta = [0 0];
        field(1).bndry(2).eta = {@(x,y) 0 , @(x,y) 0};
        % top edge
        field(1).bndry(3).loc = @(x,y) y;
        field(1).bndry(3).alpha = [0 1];
        field(1).bndry(3).beta = [1 0];
        field(1).bndry(3).eta = {@(x,y) 0 , @(x,y) 0};
        % left edge
        field(1).bndry(4).loc = @(x,y) x+L;
        field(1).bndry(4).alpha = [0 1];
        field(1).bndry(4).beta = [1 0];
        field(1).bndry(4).eta = {@(x,y) 0 , @(x,y) 0};
        
        % u_m
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
        field(3).bndry(4).loc = @(x,y) x+L;
        field(3).bndry(4).alpha = [0 0];
        field(3).bndry(4).beta = [1 1];
        field(3).bndry(4).eta = {@(x,y) 0 , @(x,y) 0};
end

%% NEWTONSOLVER------------------------------------------------------------
relTol = 1e-10;
absTol = 1e-10;
maxIter = 1;

%% NULLSPACE---------------------------------------------------------------
switch SETBC
    case 1
        field(2).nullSpace = 'true';
end

%% STREAMFUNCTION
switch SETPHI
    case 1
        existStreamFun = 'true';
end