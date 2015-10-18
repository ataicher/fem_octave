function [A, b, C] = discPhiTest(phi_m,phi_p,L)

Theta = 0;

R_m = (((3+phi_m -4*phi_m^2)/3)*phi_m^(1+2*Theta))^(-1/2);
R_p = (((3+phi_p -4*phi_p^2)/3)*phi_p^(1+2*Theta))^(-1/2);
E_m = exp(-R_m*L);
E_p = exp(-R_p*L);
S_m = (1-phi_m)*((1-4*phi_m)/(3+phi_m-4*phi_m^2))*(phi_m/R_m);
S_p = (1-phi_p)*((1-4*phi_p)/(3+phi_p-4*phi_p^2))*(phi_p/R_p);
F_m = phi_m^(2+2*Theta)*(1-phi_m);
F_p = phi_p^(2+2*Theta)*(1-phi_p);
    
    A = [     E_p                1                  0                0        0  0 ; ...
        0                 0                  1               E_m       0  0 ; ...
        1/R_p           -E_p/R_p               0                0        1  0 ; ...
        0                 0               E_m/R_m          -1/R_m      0  1 ; ...
        F_p             F_p*E_p           -F_m*E_m           -F_m       0  0 ; ...
        -phi_m*F_p*R_p  phi_m*F_p*R_p*E_p  phi_p*F_m*R_m*E_m  -phi_p*F_m*R_m  0  0 ];
    
    b = [1 1 0 0 F_p-F_m 0]';
    
    C = A\b;