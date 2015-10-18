phi_0 = .024;
w_0 = phi_0^2*4.9e-6;
U_0 = 3.18e-9;
c = w_0/U_0

x = 2*(0:.05:1);
y = 2*(0:.05:1);
[X,Y] = meshgrid(x,y);

alpha = 0;
alpha = alpha*pi/180;
A = 2*sin(alpha)^2/(pi-2*alpha-sin(2*alpha));
B = 2/(pi-2*alpha-sin(2*alpha));
r = sqrt(X.^2+Y.^2);
theta = atan(Y./X);

psi_m = r.*(A*sin(theta) - B*theta.*cos(theta))
psi_f = -c*(2*B./r + r).*sin(theta) + psi_m;

p = (-2*B./r + r).*cos(theta);

ux = (4/pi)*c*cos(2*theta)./(r.^2) - c - (2/pi)*cos(theta).^2./r + (1\pi)*theta.*sin(2*theta)./r;
uy = (4/pi)*c*sin(2*theta) - (1/pi)*sin(2*theta)./r + (2\pi)*theta.*sin(theta).^2./r;


% figure(1), surf(XX,YY,ux)
% xlabel('x'), ylabel('y')
% figure(2), surf(XX,YY,uy)
% xlabel('x'), ylabel('y')
% figure(3), quiver(XX,YY,ux,uy)
% xlabel('x'), ylabel('y')
% figure(4), surf(XX,YY,p)
% xlabel('x'), ylabel('y')
figure(5), 
hold off
v = psi_f(:,end)
contour(X,Y,psi_f,v,'--','LineWidth',2)
hold on
v = psi_m(:,end)
contour(X,Y,psi_m,v,'LineWidth',2)
legend('matrix','melt')
xlabel('depth','FontSize',20), ylabel('x','FontSize',20)
title(['melt and matrix flow lines for constant porosity phi_0 = ', num2str(phi_0),', ridge spreading rate U_0 = ', num2str(U_0)],'FontSize',20)  