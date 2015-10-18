% STOKESMANUFACTUREDSOLUTION
%
% Solve the PDE [-1,1]^2
%   -div gradU + gradP = f
%   divU = g
% for all v in V
%    (gradU , gradV) - (p, divV) - (f , v) - <n'*(gradU-pI)*n , v*n> - <tau'*gradU*n , v*tau> = 0
% for all q in Q
%    -(divU , q) + (g ,q ) = 0

%% MESH--------------------------------------------------------------------
% uniform rectangular mesh on [-1 1]^2
[nodes, edges, cells] = uniformQuadMesh(2,2,[-1 1],[-1 1]);

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

%% MANUFACTURED SOLUTION---------------------------------------------------
field(1).uExact{1} = @(x,y) x;
field(1).uExact{2} = @(x,y) y;
field(2).uExact{1} = @(x,y) 1;

field = symbolicDerivatives(field,2);

Ux = field(1).uExact{1}; Uy = field(1).uExact{2};
dxUx = field(1).gradUExact{1,1}; dyUx = field(1).gradUExact{1,2}; dxUy = field(1).gradUExact{2,1}; dyUy = field(1).gradUExact{2,2};
dxxUx = field(1).grad2UExact{1,1,1}; dyyUx = field(1).grad2UExact{1,2,2}; dxxUy = field(1).grad2UExact{2,1,1}; dyyUy = field(1).grad2UExact{2,2,2};
P = field(2).uExact{1}; 
dxP = field(2).gradUExact{1,1}; dyP = field(2).gradUExact{1,2};

%% VARIATIONAL FORM--------------------------------------------------------
fx = @(x,y) -(dxxUx(x,y) + dyyUx(x,y)) + dxP(x,y);
fy = @(x,y) -(dxxUy(x,y) + dyyUy(x,y)) + dyP(x,y);
g = @(x,y) dxUx(x,y) + dyUy(x,y);
% -(f , v)
field(1).varForm.v{1} = @(u, gradU, x, y) -fx(x,y);
field(1).varForm.v{2} = @(u, gradU, x, y) -fy(x,y);
% (gradU , gradV) - (p , divV)
field(1).varForm.gradV{1,1} = @(u, gradU, x, y) gradU.field{1}(1,1) - u.field{2}(1);
field(1).varForm.gradV{1,2} = @(u, gradU, x, y) gradU.field{1}(1,2);
field(1).varForm.gradV{2,1} = @(u, gradU, x, y) gradU.field{1}(2,1);
field(1).varForm.gradV{2,2} = @(u, gradU, x, y) gradU.field{1}(2,2) - u.field{2}(1);
% -(divU , v) - ( g , v)
field(2).varForm.v{1} = @(u, gradU, x, y) -gradU.field{1}(1,1) - gradU.field{1}(2,2) + g(x,y);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0.1;

% (gradW1, gradV1)
JacForm.field{1,1}.gradW{1,1}.gradV{1,1} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.gradW{1,2}.gradV{1,2} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.gradW{2,1}.gradV{2,1} = @(u,gradU,x,y) 1;
JacForm.field{1,1}.gradW{2,2}.gradV{2,2} = @(u,gradU,x,y) 1;
% -(w2,divV1)
JacForm.field{2,1}.w{1}.gradV{1,1} = @(u,gradU,x,y) -1;
JacForm.field{2,1}.w{1}.gradV{2,2} = @(u,gradU,x,y) -1;
% -(divW1,v2)
JacForm.field{1,2}.gradW{1,1}.v{1} = @(u,gradU,x,y) -1;
JacForm.field{1,2}.gradW{2,2}.v{1} = @(u,gradU,x,y) -1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
% MANUFACTURED SOLUTION
% bottom edge
n = [0 -1];
tau = [1 0];
alpha = [0 0];
beta = [1 1];
field(1).bndry(1).loc = @(x,y) y+1;
field(1).bndry(1).alpha = alpha;
field(1).bndry(1).beta = beta;
field(1).bndry(1).eta = {@(x,y) alpha(1)*(n(1)*dxUx(x,y)*n(1) + n(1)*dyUx(x,y)*n(2) + n(2)*dxUy(x,y)*n(1) + n(2)*dyUy(x,y)*n(2) - P(x,y)) + ...
    beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2)), ... 
    @(x,y) alpha(2)*(tau(1)*dxUx(x,y)*n(1) + tau(1)*dyUx(x,y)*n(2) + tau(2)*dxUy(x,y)*n(1) + tau(2)*dyUy(x,y)*n(2)) + ...
    beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2))};

% right edge
n = [1 0];
tau = [0 1];
alpha = [0 0];
beta = [1 1];
field(1).bndry(2).loc = @(x,y) x-1;
field(1).bndry(2).alpha = alpha;
field(1).bndry(2).beta = beta;
field(1).bndry(2).eta = {@(x,y) alpha(1)*(n(1)*dxUx(x,y)*n(1) + n(1)*dyUx(x,y)*n(2) + n(2)*dxUy(x,y)*n(1) + n(2)*dyUy(x,y)*n(2) - P(x,y)) + ...
    beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2)), ... 
    @(x,y) alpha(2)*(tau(1)*dxUx(x,y)*n(1) + tau(1)*dyUx(x,y)*n(2) + tau(2)*dxUy(x,y)*n(1) + tau(2)*dyUy(x,y)*n(2)) + ...
    beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2))};

% top edge
n = [0 1];
tau = [-1 0];
alpha = [0 0];
beta = [1 1];
field(1).bndry(3).loc = @(x,y) y-1;
field(1).bndry(3).alpha = alpha;
field(1).bndry(3).beta = beta;
field(1).bndry(3).eta = {@(x,y) alpha(1)*(n(1)*dxUx(x,y)*n(1) + n(1)*dyUx(x,y)*n(2) + n(2)*dxUy(x,y)*n(1) + n(2)*dyUy(x,y)*n(2) - P(x,y)) + ...
    beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2)), ... 
    @(x,y) alpha(2)*(tau(1)*dxUx(x,y)*n(1) + tau(1)*dyUx(x,y)*n(2) + tau(2)*dxUy(x,y)*n(1) + tau(2)*dyUy(x,y)*n(2)) + ...
    beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2))};

% left edge
n = [-1 0];
tau = [0 -1];
alpha = [0 0];
beta = [1 1];
field(1).bndry(4).loc = @(x,y) x+1;
field(1).bndry(4).alpha = alpha;
field(1).bndry(4).beta = beta;
field(1).bndry(4).eta = {@(x,y) alpha(1)*(n(1)*dxUx(x,y)*n(1) + n(1)*dyUx(x,y)*n(2) + n(2)*dxUy(x,y)*n(1) + n(2)*dyUy(x,y)*n(2) - P(x,y)) + ...
    beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2)), ... 
    @(x,y) alpha(2)*(tau(1)*dxUx(x,y)*n(1) + tau(1)*dyUx(x,y)*n(2) + tau(2)*dxUy(x,y)*n(1) + tau(2)*dyUy(x,y)*n(2)) + ...
    beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2))};

%% SOLVER
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 3;

%% NULLSPACE
% pick fields that have a constant nullSpace
field(2).nullSpace = 'true';

%% DISPLAY AND DEBUG OPTIONS
displayMesh = 'true';
% displayInitialGuess = 'true';
outputResidualError = 'true';
% displayResidualError = 'true';
displayCurrentNewtonIterate = 'true';
displayComputedSolution = 'true';
displayExactSolution = 'true';
% displayStreamLines = 'true';
% displayJacobianSpyPlot = 'true';