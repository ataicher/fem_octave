% StokesManufacturedSolutionNL.m
%
% Solve the PDE on [-1,1]^2
%   -div((1+mu*|u|^2)*gradU) + gradP = f
%   divU = g
%
% for all v in V
%    ((1+mu*|u|^2)*gradU , gradV) - (p , divV) - (f , v) - <(1+mu*|u|^2)*n'*(gradU-pI)*n , v*n> - <(1+mu*|u|^2)*tau'*gradU*n , v*tau> = 0
% for all q in Q
%    -(divU , q) + (g , q) = 0

%% MESH--------------------------------------------------------------------
% uniform rectangular mesh on [-1 1]^2
[nodes, edges, cells] = uniformQuadMesh(20,20,[-1 1],[-1 1]);
% [nodes, edges, cells] = QuadMesh([-1 -.9 -.8 -.4  1],[-1 .4 .8 .9 1]);

%% VARIABLES---------------------------------------------------------------
name.field{1} = 'velocity';
numComp.field{1} = 2;
name.field{2} = 'pressure';
numComp.field{2} = 1;

%% GUASS QUADRATURE--------------------------------------------------------
prec = 5;

%% FINITE ELEMENT BASIS----------------------------------------------------
% basisType.field{1} = @P1Quad2;
% basisType.field{2} = @piecewiseConst;
% basisType.field{1} = @AWVelocity;
% basisType.field{2} = @piecewiseConst;
basisType.field{1} = @TaylorHoodVelocity;
basisType.field{2} = @P1Quad;

%% MANUFACTURED SOLUTION
uE.field{1}{1} = @(x,y) x*y;
uE.field{1}{2} = @(x,y) x*y;
uE.field{2}{1} = @(x,y) x;
[gradUE, grad2UE] = symDerivatives(uE);
Ux = uE.field{1}{1}; Uy = uE.field{1}{2}; 
dxUx = gradUE.field{1}{1,1}; dyUx = gradUE.field{1}{1,2}; dxUy = gradUE.field{1}{2,1}; dyUy = gradUE.field{1}{2,2};
dxxUx = grad2UE.field{1}{1,1,1}; dyyUx = grad2UE.field{1}{1,2,2}; dxxUy = grad2UE.field{1}{2,1,1}; dyyUy = grad2UE.field{1}{2,2,2};
P = uE.field{2}{1}; 
dxP = gradUE.field{2}{1,1}; dyP = gradUE.field{2}{1,2};

%% VARIATIONAL FORM--------------------------------------------------------
% nonlinearity
mu = 10;
% manufactured solution forcing functions
fx = @(x,y) -2*mu*((Ux(x,y)*dxUx(x,y) + Uy(x,y)*dxUy(x,y))*dxUx(x,y) + ...
                   (Ux(x,y)*dyUx(x,y) + Uy(x,y)*dyUy(x,y))*dyUx(x,y)) - ...
            (1 + mu*(Ux(x,y)^2 + Uy(x,y)^2))*(dxxUx(x,y) + dyyUx(x,y)) + ...
            dxP(x,y);
fy = @(x,y) -2*mu*((Ux(x,y)*dxUx(x,y) + Uy(x,y)*dxUy(x,y))*dxUy(x,y) + ...
                   (Ux(x,y)*dyUx(x,y) + Uy(x,y)*dyUy(x,y))*dyUy(x,y)) - ...
            (1 + mu*(Ux(x,y)^2 + Uy(x,y)^2))*(dxxUy(x,y) + dyyUy(x,y)) + ...
            dyP(x,y);
g = @(x,y) -dxUx(x,y) - dyUy(x,y);
% (f , v)
varForm.field{1}.v{1} = @(u, gradU, x, y) -fx(x,y);
varForm.field{1}.v{2} = @(u, gradU, x, y) -fy(x,y);
% ((1+mu*|u|^2)*gradU , gradV) - (p , divV) 
varForm.field{1}.gradV{1,1} = @(u, gradU, x, y) (1 + mu*(u.field{1}(1)^2 + u.field{1}(2)^2))*gradU.field{1}(1,1) - u.field{2};
varForm.field{1}.gradV{1,2} = @(u, gradU, x, y) (1 + mu*(u.field{1}(1)^2 + u.field{1}(2)^2))*gradU.field{1}(1,2);
varForm.field{1}.gradV{2,1} = @(u, gradU, x, y) (1 + mu*(u.field{1}(1)^2 + u.field{1}(2)^2))*gradU.field{1}(2,1);
varForm.field{1}.gradV{2,2} = @(u, gradU, x, y) (1 + mu*(u.field{1}(1)^2 + u.field{1}(2)^2))*gradU.field{1}(2,2) - u.field{2};
% -(divU, q)
varForm.field{2}.v{1} = @(u, gradU, x, y) -gradU.field{1}(1,1) - gradU.field{1}(2,2) - g(x,y);

%% JACOBAIN FORM-----------------------------------------------------------
% finite difference parameter
h = 0;
% ((1+mu*|u|^2)*gradW , gradV)
JacForm.field{1,1}.gradW{1,1}.gradV{1,1} = @(u,gradU,x,y) 1 + mu*(u.field{1}(1)^2 + u.field{1}(2)^2);
JacForm.field{1,1}.gradW{1,2}.gradV{1,2} = @(u,gradU,x,y) 1 + mu*(u.field{1}(1)^2 + u.field{1}(2)^2);
JacForm.field{1,1}.gradW{2,1}.gradV{2,1} = @(u,gradU,x,y) 1 + mu*(u.field{1}(1)^2 + u.field{1}(2)^2);
JacForm.field{1,1}.gradW{2,2}.gradV{2,2} = @(u,gradU,x,y) 1 + mu*(u.field{1}(1)^2 + u.field{1}(2)^2);
% ((2*mu*(ux+uy)*gradU)*w , gradV)
JacForm.field{1,1}.w{1}.gradV{1,1} = @(u,gradU,x,y) 2*mu*u.field{1}(1)*gradU.field{1}(1,1);
JacForm.field{1,1}.w{2}.gradV{1,1} = @(u,gradU,x,y) 2*mu*u.field{1}(2)*gradU.field{1}(1,1);
JacForm.field{1,1}.w{1}.gradV{1,2} = @(u,gradU,x,y) 2*mu*u.field{1}(1)*gradU.field{1}(1,2);
JacForm.field{1,1}.w{2}.gradV{1,2} = @(u,gradU,x,y) 2*mu*u.field{1}(2)*gradU.field{1}(1,2);
JacForm.field{1,1}.w{1}.gradV{2,1} = @(u,gradU,x,y) 2*mu*u.field{1}(1)*gradU.field{1}(2,1);
JacForm.field{1,1}.w{2}.gradV{2,1} = @(u,gradU,x,y) 2*mu*u.field{1}(2)*gradU.field{1}(2,1);
JacForm.field{1,1}.w{1}.gradV{2,2} = @(u,gradU,x,y) 2*mu*u.field{1}(1)*gradU.field{1}(2,2);
JacForm.field{1,1}.w{2}.gradV{2,2} = @(u,gradU,x,y) 2*mu*u.field{1}(2)*gradU.field{1}(2,2);
% -(w,divV)
JacForm.field{2,1}.w{1}.gradV{1,1} = @(u,gradU,x,y) -1;
JacForm.field{2,1}.w{1}.gradV{2,2} = @(u,gradU,x,y) -1;
% -(divW,v)
JacForm.field{1,2}.gradW{1,1}.v{1} = @(u,gradU,x,y) -1;
JacForm.field{1,2}.gradW{2,2}.v{1} = @(u,gradU,x,y) -1;

%% BOUNDARY CONDITIONS-----------------------------------------------------
% MANUFACTURED SOLUTION
% bottom edge
loc = @(x,y) y+1;
n = [0 -1];
tau = [1 0];
alpha = [0 0];
beta = [1 1];
eta{1} = @(x,y) alpha(1)*(1 + mu*(Ux(x,y)^2 + Uy(x,y)^2))* ...
    (n(1)*dxUx(x,y)*n(1) + n(1)*dyUx(x,y)*n(2) + n(2)*dxUy(x,y)*n(1) + n(2)*dyUy(x,y)*n(2) - P(x,y)) + ...
    beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2));
eta{2} = @(x,y) alpha(2)*(1 + mu*(Ux(x,y)^2 + Uy(x,y)^2))* ...
    (tau(1)*dxUx(x,y)*n(1) + tau(1)*dyUx(x,y)*n(2) + tau(2)*dxUy(x,y)*n(1) + tau(2)*dyUy(x,y)*n(2)) + ...
    beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2));
BC.field{1}.bndry{1} = assembleBoundary(loc,alpha,beta,eta);
% right edge
loc = @(x,y) x-1;
n = [1 0];
tau = [0 1];
alpha = [0 0];
beta = [1 1];
eta{1} = @(x,y) alpha(1)*(1 + mu*(Ux(x,y)^2 + Uy(x,y)^2))* ...
    (n(1)*dxUx(x,y)*n(1) + n(1)*dyUx(x,y)*n(2) + n(2)*dxUy(x,y)*n(1) + n(2)*dyUy(x,y)*n(2) - P(x,y)) + ...
    beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2));
eta{2} = @(x,y) alpha(2)*(1 + mu*(Ux(x,y)^2 + Uy(x,y)^2))* ...
    (tau(1)*dxUx(x,y)*n(1) + tau(1)*dyUx(x,y)*n(2) + tau(2)*dxUy(x,y)*n(1) + tau(2)*dyUy(x,y)*n(2)) + ...
    beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2));
BC.field{1}.bndry{2} = assembleBoundary(loc,alpha,beta,eta);
% top edge
loc = @(x,y) y-1;
n = [0 1];
tau = [-1 0];
alpha = [0 0];
beta = [1 1];
eta{1} = @(x,y) alpha(1)*(1 + mu*(Ux(x,y)^2 + Uy(x,y)^2))* ...
    (n(1)*dxUx(x,y)*n(1) + n(1)*dyUx(x,y)*n(2) + n(2)*dxUy(x,y)*n(1) + n(2)*dyUy(x,y)*n(2) - P(x,y)) + ...
    beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2));
eta{2} = @(x,y) alpha(2)*(1 + mu*(Ux(x,y)^2 + Uy(x,y)^2))* ...
    (tau(1)*dxUx(x,y)*n(1) + tau(1)*dyUx(x,y)*n(2) + tau(2)*dxUy(x,y)*n(1) + tau(2)*dyUy(x,y)*n(2)) + ...
    beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2));
BC.field{1}.bndry{3} = assembleBoundary(loc,alpha,beta,eta);
% left edge
loc = @(x,y) x+1;
n = [-1 0];
tau = [0 -1];
alpha = [0 0];
beta = [1 1];
eta{1} = @(x,y) alpha(1)*(1 + mu*(Ux(x,y)^2 + Uy(x,y)^2))* ...
    (n(1)*dxUx(x,y)*n(1) + n(1)*dyUx(x,y)*n(2) + n(2)*dxUy(x,y)*n(1) + n(2)*dyUy(x,y)*n(2) - P(x,y)) + ...
    beta(1)*(Ux(x,y)*n(1) + Uy(x,y)*n(2));
eta{2} = @(x,y) alpha(2)*(1 + mu*(Ux(x,y)^2 + Uy(x,y)^2))* ...
    (tau(1)*dxUx(x,y)*n(1) + tau(1)*dyUx(x,y)*n(2) + tau(2)*dxUy(x,y)*n(1) + tau(2)*dyUy(x,y)*n(2)) + ...
    beta(2)*(Ux(x,y)*tau(1) + Uy(x,y)*tau(2));
BC.field{1}.bndry{4} = assembleBoundary(loc,alpha,beta,eta);

%% SOLVER------------------------------------------------------------------
% Solver uses Newton iteration
relTol = 1e-10;
absTol = 1e-10;
maxIter = 10;

%% NULLSPACE---------------------------------------------------------------
% pick fields that have a constant nullSpace
nullSpace.field{2} = 'true';

%% DISPLAY AND DEBUG OPTIONS-----------------------------------------------
% displayMesh = 'true';
% displayInitialGuess = 'true';
% displayResidualError = 'true';
% displayCurrentNewtonIterate = 'true';
displayComputedSolution = 'true';
exactSolutionExist = 'true';
displayExactSolution = 'true';
% displayStreamLines = 'true';
% displayJacobianSpyPlot = 'true';