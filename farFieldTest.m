% sigma11 = @(x,y) 2*dxUmx(x,y)-2/3*(dxUmx(x,y)+dyUmy(x,y))-Qm(x,y);
% sigma12 = @(x,y) dyUmx(x,y)+dxUmy(x,y);
% sigma21 = @(x,y) dyUmx(x,y)+dxUmy(x,y);
% sigma22 = @(x,y) 2*dyUmy(x,y)-2/3*(dxUmx(x,y)+dyUmy(x,y))-Qm(x,y);
% sigma = @(x,y) [ sigma11(x.y) sigma12(x,y) ; sigma21(x,y) sigma22(x,y)];
% t=0:.01:Lx;
% % sigma = @(x,y) dxUmx(x,y) - 2/3*(dxUmx(x,y)-dyUmy(x,y));
% sigma = @(x,y) [ 2*dxUmx(x,y)-2/3*(dxUmx(x,y)-dyUmy(x,y))-Qm(x,y)  -dyUmx(x,y)-dxUmy(x,y) ; ...
%     -dyUmx(x,y)-dxUmy(x,y) 2*dyUmy(x,y)-2/3*(dxUmx(x,y)-dyUmy(x,y))-Qm(x,y)];
% sigmaVec = zeros(size(t));
% for i=1:length(t)
%     A = sigma(t(i),Ly);
%     n = [0 1]';
%     tan = [-1 0]';
% sigmaVecN(i) = n'*A*n + t(i);
% sigmaVecTan(i) = tan'*A*n;
% end
% figure
% plot(t,sigmaVecN)
% figure
% plot(t,sigmaVecTan)