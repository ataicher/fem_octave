L = 1;
Theta = 0;
phi_m = .01;
phi_p = .04;
R_m = ((1-phi_m)*((4/3)*phi_m+1)*phi_m^(1+2*Theta))^(-1/2);
R_p = ((1-phi_p)*((4/3)*phi_p+1)*phi_p^(1+2*Theta))^(-1/2);
E_m = exp(-R_m*L);
E_p = exp(-R_p*L);
phiFM_p = phi_p^(2*(1+Theta))*(1-phi_p);
phiFM_m = phi_m^(2*(1+Theta))*(1-phi_m);
phiFP_p = phi_p^(1+2*Theta)*(1-phi_p)^2;
phiFP_m = phi_m^(1+2*Theta)*(1-phi_m)^2;

A = [E_p         1            0          0    0 0 ; ...
      0          0            1         E_m   0 0 ; ...
     1/R_p    -E_p/R_p        0          0    1 0 ; ...
      0          0         E_m/R_m     -1/R_m 0 1 ; ...
   phiFM_p  phiFM_p*E_p -phiFM_m*E_m -phiFM_m 0 0 ; ...
   -phiFP_p*R_p+(1-phi_p)/R_p      phiFP_p*R_p*E_p-(1-phi_p)*E_p/R_p     phiFP_m*R_m*E_m-(1-phi_m)*E_m/R_m     -phiFP_m*R_m+(1-phi_m)/R_m    1-phi_p    -(1-phi_m)];

b = [1 1 0 0 phiFM_p-phiFM_m 0]';


C = A\b;

C1_p = C(1);
C2_p = C(2);
C1_m = C(3);
C2_m = C(4);
CPInt_p = C(5);
CPInt_m = C(6);

x_m = -L:.001:0;
x_p = 0:.001:L;


vr_m = phi_m^(1+Theta)*(1-phi_m)*(1 - C1_m*exp(-R_m*(x_m+L)) - C2_m*exp(R_m*x_m));
vr_p = phi_p^(1+Theta)*(1-phi_p)*(1 - C1_p*exp(-R_p*x_p) - C2_p*exp(R_p*(x_p-L)));
p_m = (1-phi_m)*(x_m + C1_m*exp(-R_m*(x_m+L))/R_m - C2_m*exp(R_m*x_m)/R_m + CPInt_m);
p_p = (1-phi_p)*(x_p + C1_p*exp(-R_p*x_p)/R_p - C2_p*exp(R_p*(x_p-L))/R_p + CPInt_p);
vm_m = phi_m^(2*(1+Theta))*(1-phi_m)*(1 - C1_m*exp(-R_m*(x_m+L)) - C2_m*exp(R_m*x_m));
vm_p = phi_p^(2*(1+Theta))*(1-phi_p)*(1 - C1_p*exp(-R_p*x_p) - C2_p*exp(R_p*(x_p-L)));
phiV_m = ones(size(x_m))*phi_m;
phiV_p = ones(size(x_p))*phi_p;
subplot(1,4,1)
plot([x_m x_p],[vr_m vr_p])
view(90,90)
title('relative velocity')
subplot(1,4,2)
plot([x_m x_p],[p_m p_p])
view(90,90)
title('fluid pressure')
subplot(1,4,3)
plot([x_m x_p],[vm_m vm_p])
view(90,90)
title('matrix velocity')
subplot(1,4,4)
plot([x_m x_p],[phiV_m phiV_p])
axis([ -L L 0 phi_p*1.1])
view(90,90)
title('porosity')