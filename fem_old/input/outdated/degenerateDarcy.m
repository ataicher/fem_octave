% Solves a second order elliptic equation with a degenerate porosity using
% lowest order RT finite elements

% EQUATIONS
% Given the porosity phi (0 <= phi < 1), find (u,p) such that
%   u = - phi \grad p
%   visc \div(phi u) + phi p = phi source
%   phi^{1/2} p = bcD    on the boundary
% Energy estimates show stability for u and phi^{1/2} p.  We believe that
% p is a poor variable when phi=0, so we define
%   q = phi^{1/2} p
% And reformulate as
%   u = - phi \grad (\phi^{-1/2} q) = 0
%   phi^{-1/2}\div(phi u) + q/visc = phi^{1/2} source/visc
%   q = bcD    on the boundary
% u       is the (unknown) velocity
% q       is the (unknown) scaled pressure
% source  is the external source (assumed constant per grid element)
% bcD     is the external boundary scaled pressure

% VARIATIONAL FORM
% Let V include DOFs on the boundary.
% Find (u,p) in V x W such that
%   a(u,v) - b(q,v) = -<bcD,phi^{(n-1)/2} v.nu>  for all v <--> u
%   b(w,u) + c(q,w) = (phi^{1/2} source,w)/visc  for all w <--> q
% where
% a(u,v) = (u,v)
% b(w,v) = (w, phi^{-1/2}\div(phi^{n/2} v) )
% c(q,w) = (q,w)/visc

% u       is the (unknown) velocity
% q       is the (unknown) scaled pressure
% source  is the external source (assumed constant per grid element)
% bcD     is the external boundary scaled pressure

% MIXED FE METHOD
% We use a rectangular domain and a rectangular grid, and lowest order RT_0
% finite elements V = Q_{1,0} X Q_{0,1}, W = Q_{0,0}.
% 
% We approximate phi as an RT0 potential in Q_{2,0} + Q_{0,2}.  It has
% edge (e) DOFs, being the average on the edge, and the element (E) average,
% denoted as phi_e and phi_c, respectively.

% Thus
%   a(v_i,v_j) = (v_i,v_j) = (1/2)( |E_i| + |E_j| ) \delta_{ij}
%   b(w_j,v_i) = 0    if phi_c = 0
%              = \sum_{e of E_j} phi_0^{-1/2} phi_e |e|  if phi_c != 0
%   c(w_i,w_j) = (w_i,w_j)/visc = |E_i|/visc \delta_{ij}

% V x W = RT_0

% //
% //  LINEAR SYSTEM
% //     /  A    -B \  /  u  \   / a \.
% //     \  B^T   C /  \  q  / = \ b /
% //    where A and C are diagonal, and
% //      a_j = -<bcD,phi^{1/2} v_j.nu> = (+/-)phi_e^{1/2} bcD |e| on boundary (or 0)
% //      b_j = (phi^{1/2} source,w_j)/visc = phi_0^{1/2} source/visc |E_j|
% //    The Shur complement for q is
% //      u = A^{-1} (B q + a)

%% PROBELM PARAMETERS
Theta = 0;
visc = 1;

% mesh
SETMESH = 1;
% rectangular mesh on [-1 1]x[-1 1]
% 1. uniform mesh with nx x ny elements
if ~exist('nx','var') || ~exist('ny','var')
    nx = 10;
    ny = 20;
end

% porosity
if ~exist('SETPHI','var')
    SETPHI = 6;
end
%   1. constant: 
%       phi = phi_0 
%   2. Polynomial equal to zero on boundaries
%       phi = (1-x^2)*(1-y^2)
%   3. Polynomial
%       phi = x^2*y^2
%   4. Non-zero only for x,y positive and goes to zero like sqrt and square
%       phi = { sqrt(x)*y^2  x,y >= 0
%             { 0             else
%   5. Quadratic in x-direction only and vanishes at x = 0
%       phi = { x^alpha  x >= 0
%             { 0        else
%   6. 
%       phi = { (x+3/4)^alpha*(y+3/4)^(2*alpha)  x,y >= -3/4
%       phi = { 0                                   else
phi_0 = .0001;
alpha = .25;

% pressure
% pressure for manufactured solution. Set SETPRES to 0 when not considering
% a manufactured solution
if ~exist('SETPRES','var')
SETPRES = 6;
end
%   1. Mix of polynomial and trig function
%       p = x^2*cos(y)
%   2. Polynomial that is zero on the boundaries
%       p = (1-x^2)^2*(1-y^2)^2
%   3. Highly varying trig function
%       p = cos(6*x^2*y)
%   4. Exact solution when SETPHI = 4 and alpha = 2, i.e.,
%       phi = { x^2  x >= 0
%             { 0     else
%      for source = x^beta
%       p = (beta*x^r1-r1*x^beta)/(r1*(beta-r1)*(beta-r2))
%   5. Mix of polynomial and trig function
%       p = x^2*cos(y)+y^2*cos(x)
%   6. Highly varying trig function
%       p = cos(6*x*y^2)
beta = 2;

% source
% source term only applied when SETPRES = 0
SETSOURCE = 0;

% under-integration
SETUNDERINTEGRATE = 1;

%% MESH--------------------------------------------------------------------
% uniform rectangular mesh on [-1 1]x[-1 1]
switch SETMESH
    
    % UNIFORM
    case 1
        Lx = [-1 1];
        Ly = [-1 1];
end

%% VARIABLES---------------------------------------------------------------
field(1).name = 'velocity';
field(1).numComp = 2;
field(2).name = 'scaled pressure';
field(2).numComp = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 3;

%% FINITE ELEMENT BASIS----------------------------------------------------
field(1).basisType = @RT0Velocity;
field(2).basisType = @piecewiseConst;

%% DEGENERATE--------------------------------------------------------------
degenerate = 'true';
if SETUNDERINTEGRATE
    underIntegrate = 'true';
end
RT0VelocityField = 1;
RT0PressureField = 2;

switch SETPHI
    case 1
        phi = @(x,y) .04;
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
        phi = @(x,y) x^alpha*(x>=0);
        dxPhi = @(x,y) alpha*x^(alpha-1)*(x>=0);
        dyPhi = @(x,y) 0;
    case 6
        phi = @(x,y) ((x+3/4)*(y+3/4)^2)^alpha*(x>=-3/4)*(y>=-3/4);
        dxPhi = @(x,y) alpha*((x+3/4)*(y+3/4)^2)^(alpha-1)*(y+3/4)^2*(x>=-3/4)*(y>=-3/4);
        dyPhi = @(x,y) alpha*((x+3/4)*(y+3/4)^2)^(alpha-1)*(x+3/4)*2*(y+3/4)*(x>=-3/4)*(y>=-3/4);
end

%% EXACT SOLUTION----------------------------------------------------------
switch SETPRES
    case 1
        p = @(x,y) x^2*cos(y);
        dxP = @(x,y) 2*x*cos(y);
        dxxP = @(x,y) 2*cos(y);
        dyP = @(x,y) -x^2*sin(y);
        dyyP = @(x,y) -x^2*cos(y);
        dxyP = @(x,y) -2*x*sin(y);
    case 2
        p = @(x,y) (1-x^2)^2*(1-y^2)^2;
        dxP = @(x,y) -4*x*(1-x^2)*(1-y^2)^2;
        dxxP = @(x,y) -4*(1-3*x^2)*(1-y^2)^2;
        dyP = @(x,y) -(1-x^2)^2*4*y*(1-y^2);
        dyyP = @(x,y) -(1-x^2)^2*4*(1-3*y^2);
        dxyP = @(x,y) 16*x*y*(1-x^2)*(1-y^2);
    case 3
        p = @(x,y) cos(6*x^2*y);
        dxP = @(x,y) -12*x*y*sin(6*x^2*y);
        dxxP = @(x,y) -12*y*sin(6*x^2*y) - 144*x^2*y^2*cos(6*x^2*y);
        dyP = @(x,y) -6*x^2*sin(6*x^2*y);
        dyyP = @(x,y) -36*x^4*cos(6*x^2*y);
        dxyP = @(x,y) -12*x*cos(6*x^2*y) - 72*x^3*y*cos(6*x^2*y);
    case 4
        r1 = (-3+sqrt(13))/2; r2 = (-3-sqrt(13))/2;
        p = @(x,y) (beta*x^r1 - r1*x^beta)/(r1*(beta-r1)*(beta-r2))*(x>=0);
        dxP = @(x,y) beta*(x^(r1-1) - x^(beta-1))/((beta-r1)*(beta-r2))*(x>=0);
        dxxP = @(x,y) beta*( (r1-1)*x^(r1-2) - (beta-1)*x^(beta-2))/((beta-r1)*(beta-r2))*(x>=0);
        dyP = @(x,y) 0;
        dyyP = @(x,y) 0;
        dxyP = @(x,y) 0;
    case 5
        p = @(x,y) x^2*cos(y) + y^2*cos(x);
        dxP = @(x,y) 2*x*cos(y) - y^2*sin(x);
        dxxP = @(x,y) 2*cos(y) - y^2*cos(x);
        dyP = @(x,y) -x^2*sin(y) + 2*y*cos(x);
        dyyP = @(x,y) -x^2*cos(y) + 2*cos(x);
        dxyP = @(x,y) -2*x*sin(y) -2*y*sin(x);
    case 6
        p = @(x,y) cos(6*x*y^2);
        dxP = @(x,y) -6*y^2*sin(6*y^2*x);
        dxxP = @(x,y) -36*y^4*cos(6*y^2*x);
        dyP = @(x,y) -12*x*y*sin(6*y^2*x);
        dyyP = @(x,y) -12*x*sin(6*y^2*x) - 144*x^2*y^2*cos(6*y^2*x);
        dxyP = @(x,y) -12*y*sin(6*x*y^2) - 72*x*y^3*cos(6*x*y^2);
end

if SETPRES~=0
    Vx = @(x,y) -phi(x,y)^(1+Theta)*dxP(x,y);
    Vy = @(x,y) -phi(x,y)^(1+Theta)*dyP(x,y);
    
    dxVx = @(x,y) -(1+Theta)*phi(x,y)^Theta*dxPhi(x,y)*dxP(x,y) - phi(x,y)^(1+Theta)*dxxP(x,y);
    dyVx = @(x,y) -(1+Theta)*phi(x,y)^Theta*dyPhi(x,y)*dxP(x,y) - phi(x,y)^(1+Theta)*dxyP(x,y);
    dxVy = @(x,y) -(1+Theta)*phi(x,y)^Theta*dyPhi(x,y)*dyP(x,y) - phi(x,y)^(1+Theta)*dxyP(x,y);
    dyVy = @(x,y) -(1+Theta)*phi(x,y)^Theta*dyPhi(x,y)*dyP(x,y) - phi(x,y)^(1+Theta)*dyyP(x,y);
   
    Q = @(x,y) sqrt(phi(x,y))*p(x,y);

    source = @(x,y) - visc*phi(x,y)*(dxxP(x,y) + dyyP(x,y)) - 2*visc*(dxPhi(x,y)*dxP(x,y) + dyPhi(x,y)*dyP(x,y)) + p(x,y); 
    
    field(1).uExact{1} = @(x,y) Vx(x,y);
    field(1).gradUExact{1,1} = @(x,y) dxVx(x,y);
    field(1).gradUExact{1,2} = @(x,y) dyVx(x,y);
    field(1).uExact{2} = @(x,y) Vy(x,y);
    field(1).gradUExact{2,1} = @(x,y) dxVy(x,y);
    field(1).gradUExact{2,2} = @(x,y) dyVy(x,y);
    
    field(2).uExact{1} = @(x,y) Q(x,y);
    field(2).gradUExact{1,1} = @(x,y) 0;
    field(2).gradUExact{1,2} = @(x,y) 0;
end

if SETPRES == 0
    switch SETSOURCE
        case 1
            source = @(x,y) -y;
        case 2
            source = @(x,y) x*y;
        case 3
            source = @(x,y) abs(x)^beta;
        case 4
            source = @(x,y) x^beta*(x>=0);
        case 5
            
    end
end

%% VARIATIONAL FORM--------------------------------------------------------
% (u , v) - (q , phi^(-1/2)*div(phi^(n/2)*v) = 0
if ~SETUNDERINTEGRATE
    field(1).varForm.v{1} = @(u, gradU, x, y) u.field{1}(1);
    field(1).varForm.v{2} = @(u, gradU, x, y) u.field{1}(2);
end

% (phi^(-1/2)*div(phi^(n/2)*u , w) + (q , w) - (phi^(1/2)*source , w) = 0
field(2).varForm.v{1} = @(u, gradU, x, y) u.field{2}(1) - phi(x,y)^(1/2)*source(x,y)/visc;

%% JACOBAIN FORM-----------------------------------------------------------
h = 0.1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
field(1).bndry(1).loc = @(x,y) (x-1)*(y-1)*(x+1)*(y+1);
field(1).bndry(1).alpha = [1 1];
field(1).bndry(1).beta = [0 0];
field(1).bndry(1).eta = {@(x,y) phi(x,y)^(1+Theta)*p(x,y), @(x,y) 0};

%% SOLVER------------------------------------------------------------------
relTol = 1e-10;
absTol = 1e-10;
maxIter = 1;

%% NULLSPACE---------------------------------------------------------------
% field(2).nullSpace = 'true';