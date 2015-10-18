% compactingColumn.m
%
% 1D compacting column, in non-dimensionalized form, on [-L 0].  The
% multidimensional mechanical equations in this case reduce to the 1D
% equations:
%
% 1. u = -phi^2(1+Theta)*q_f'
% 2. u' = -phi(q_f-q_m)
% 3. (q_m - 1/3*(1-4*phi)*v_m')' = 1-phi
% 4. v_m' = phi*(q_f-q_m) = 0
%
% with boundary conditions:
% 1.
%              u = 0
%             v_m = 0
%               ___
%                |
%                |
%                |
%                |
%               _|_
%
%              u = 0
%             v_m = 0
%
%
% 2.
%             q_f = 0
%             v_m = 0
%               ___
%                |
%                |
%                |
%                |
%               _|_
%
%             u = 0
%             v_m = 0
%
% To replicate this problem in 2D we solve the full system on the
% rectangular domain [-L L]x[0 1]:
% Set sigma(v_m) = 2*(1-phi)*grad v_m - (5-2*phi)/3 * div v_m

% I. Standard Method
% 1. (phi^(-2-2*Theta))u , Psi) - (q_f , div Psi) + <q_f , Psi*n> = 0
%
% 2. (div(u) , w_f) + ( phi(q_f-q_m) , w_f) = 0
%
% 3. -(q_m , div Psi_m) + (div sigma(v_m) , grad Psi_m) + <sigma(v_m) , Psi_m*n> - ( (1-phi)*drho*g , Psi_m) = 0
%
% 4. (div v_m  , w_m) - ( phi*(q_f-q_m) , w) = 0
%
% II. Auxilliary Velocity
% Set phi^(2+2*Theta)*u_A = u
%
% 1. (u , Psi_A) -(phi^(2+2*Theta)*u_A , Psi_A) = 0
%
% 2. (u_A , Psi) - (q_f , div Psi) + <q_f , Psi*n> = 0
%
% 3. (div(u) , w_f) + ( phi(q_f-q_m) , w_f) = 0
%
% 4. -(q_m , div Psi_m) + (div sigma(v_m) , grad Psi_m) + <sigma(v_m) , Psi_m*n> - ( (1-phi)*drho*g , Psi_m) = 0
%
% 5. (div v_m  , w_m) - ( phi*(q_f-q_m) , w) = 0
%
% III. Symmetric method
% Set u = phi^(1+Theta)*u
%
% 1. (u , Psi) - (q_f , div phi^(1+Theta)*Psi) + <q_f , Psi*n> = 0
%
% 2. (div(phi^(1+Theta)*u) , w_f) + ( phi(q_f-q_m) , w_f) = 0
%
% 3. -(q_m , div Psi_m) + (div sigma(v_m) , grad Psi_m) + <sigma(v_m) , Psi_m*n> - ( (1-phi)*drho*g , Psi_m) = 0
%
% 4. (div v_m  , w_m) - ( phi*(q_f-q_m) , w) = 0
%
% IV. Scaled formulation
% Set u = phi^(1+Theta)*u
%     q_f = phi^(1/2)*q_f
%
% 1. (u , Psi) - (q_f , phi^(-1/2)*div phi^(1+Theta)*Psi) + <phi^(1/2+Theta)*q_f , Psi*n> = 0
%
% 2. (phi^(-1/2)* div(phi^(1+Theta)*u) , w_f) + (q_f , w_f) - (phi^(1/2)*q_m , w_f) = 0
%
% 3. -(q_m , div Psi_m) + (div sigma(v_m) , grad Psi_m) + <sigma(v_m) , Psi_m*n> - ( (1-phi)*drho*g , Psi_m) = 0
%
% 4. (div v_m  , w_m) - ( phi*(q_f-q_m) , w) = 0
%
% with boundary condtions:
% 1.
%                       u = 0
%                      v_m = 0
%                   _____________
%                  |             |
%                  |             |
%                  |             |
%       u*n = 0    |             |    u*n = 0
%      v_m*n = 0   |             |   v_m*n = 0
%   tau'*v_m*n = 0 |             | tau'*v_m*n = 0
%                  |             |
%                  |             |
%                  |_____________|
%
%                       u = 0
%                      v_m = 0
%
% 2.
%                      q_f = 0
%                      v_m = 0
%                   _____________
%                  |             |
%                  |             |
%                  |             |
%       u*n = 0    |             |    u*n = 0
%      v_m*n = 0   |             |   v_m*n = 0
%   tau'*v_m*n = 0 |             | tau'*v_m*n = 0
%                  |             |
%                  |             |
%                  |_____________|
%
%                      u = 0
%                      v_m = 0

%% PROBELM PARAMETERS
Theta = 0;

% mesh
% rectangular mesh on [-L L]x[0 1]
% 1. uniform mesh with nx x 1 elements
% 2. mesh refined around x = 0
SETMESH = 1;
L = 1;
if ~exist('nx','var')
    nx = 100;
end

% porosity
%   1. constant: phi = phi_0
%   2. Gaussian: phi = A*exp(-x^2/B^2) + phi_0
%   3. Quadratic:
%             { phi_m               x <= 0
%       phi = { phi_m + x^2   0 < x <= (phi_p-phi_m)^(1/2)
%             { phi_p               x > (phi_p-phi_m)^(1/2)
%   4. Square root:
%             { phi_m                   x <= 0
%       phi = { phi_m + x^(1/2)   0 < x <= (phi_p-phi_m)^2
%             { phi_p                   x > ((phi_p-phi_m)/D)^2
%   5. Discontinuous jump:
%       phi = {phi_m x <  0
%             {phi_p x >= 0
if ~exist('SETPHI','var')
    SETPHI = 5;
end
phi_0 = 0.04;
if ~exist('phi_m','var')
    phi_m = 0;
end
if ~exist('phi_p','var')
    phi_p = .04;
end
A = .04;
B = .2;

% boundary condition
%   1. no flow at both top and bottom boundaries
%   2. zero fluid pressure at top boundary and no flow at bottom
SETBC = 1;

% pick implemenatation
%   1. standard (divide by phi^(2+2*Theta))
%   2. auxilliary velocity
%   3. symmetric
%   4. scaled
if ~exist('FORMULATION','var')
    FORMULATION = 4;
end

%% POROSITY
switch SETPHI
    % constant
    case 1
        phi = @(x,y) phi_0;
        dxPhi = @(x,y) 0;      
        % Gaussian
    case 2
        phi = @(x,y) A*exp(-x^2/(B^2)) + phi_0;
        dxPhi = @(x,y) -(2*A*x/(B^2))*exp(-x^2/(B^2));
        
        % quadratic
    case 3
        phi = @(x,y) phi_m*(x<0) + (phi_m + x^2)*(x>=0)*(x<=sqrt(phi_p-phi_m))  + phi_p*(x > sqrt(phi_p-phi_m));
        dxPhi = @(x,y) 2*x*(x>0)*(x<=sqrt(phi_p-phi_m));
        
        % square root
    case 4
        phi = @(x,y) phi_m*(x<0) + (phi_m + (phi_p-phi_m)*sqrt(x/L))*(x>=0);
        dxPhi = @(x,y) 1/(2*sqrt(x))*(x>0)*(x<=(phi_p-phi_m)^2);
        
        % discontinous jump at x = 0
    case 5
        phi = @(x,y) phi_m*(x<0) + phi_p*(x>=0);
        dxPhi = @(x,y) 0;
        phi = @(x,y) phi_m*(x<2/nx) + (phi_m + nx*(phi_p-phi_m)*(x+2/nx)/4 )* (x>=-2/nx)*(x<=2/nx) + phi_p*(x>2/nx);
        dxPhi = @(x,y) 3*(phi_p-phi_m)/nx;
end
dyPhi = @(x,y) 0;

%% EXACT SOLUTION
if SETBC == 1
    switch SETPHI
        % constant
        case 1
            constPhiExactSol;
            % quadratic approx
        case 3
            if L^2 - phi_p > 1e-10 || Theta ~= 0 || phi_p > .2
                warning('no approximate solution when phi > phi_p, phi > .2, or Theta ~= 0')
            else
                sqrPhiApproxSol;
            end
            % discontinuous porosity
        case 5
            discPhiExactSol;
    end
end

%% MESH--------------------------------------------------------------------
switch SETMESH
    
    % uniform mesh
    case 1
        ny = 1;
        Lx = [-L L]; Ly = [0 1];
        
        % refined at center
    case 2
        x = [-L:L/50:-L/10 -L/10+L/1000:L/1000:L/10-L/1000 L/10:L/50:L];
        y = [0 1];
end

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
    SETUNDERINTEGRATE = 0;
    if SETUNDERINTEGRATE
        underIntegrate = 'true';
    end
    RT0VelocityField = 1;
    RT0PressureField = 2;
else
    SETDEGENERATE = 0;
end

%% VARIATIONAL FORM--------------------------------------------------------
mantleMechanicsVarForm;

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
switch SETBC
    
    % no flow at top and bottom
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
        
        % pressure at top and no flow at bottom
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
if SETBC == 1
    field(2).nullSpace = 'true';
%     field(4).nullSpace = 'true';
end

%% STREAMFUNCTION
if SETPHI == 1
    existStreamFun = 'true';
end