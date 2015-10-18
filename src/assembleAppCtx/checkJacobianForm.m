function checkJacobianForm(appCtx)

globalSize = appCtx.globalSize;
h = appCtx.h;

if h
    error('myApp:argChk','Jacobian form is not defined.')
end

% random guess
u = rand(globalSize,1);

% compute exact Jacobian
appCtx.h = 0;
Jac = computeJacobian(u,appCtx);

% compute series of finite difference Jacobians
h=10.^-(0:6);
FDError = zeros(size(h));
for i=1:6
    appCtx.h = 1e-i;
    JacFD = computeFiniteDiffJacobian(u,appCtx);
    FDError(i) = norm(Jac-JacFD);    
end

loglog(1./h,FDError,'b*-')
xlabel('1/h')
ylabel('norm')
title(['convergence of FD Jacobian to user defined Jacobian form - problem size: ', num2str(globalSize), ' DOF'])
