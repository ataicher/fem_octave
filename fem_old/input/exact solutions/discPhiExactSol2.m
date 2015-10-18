% check for input variables make sense
if ~(exist('phi_m', 'var') && exist('phi_p', 'var') && exist('Theta', 'var') && exist('L', 'var'))
    error('phi_m, phi_p, Theta, or L does not exist')
end
if ~(isnumeric(phi_m) && phi_m >= 0 && phi_m < 1)
    error('phi_m must be a number in [0,1)')
end
if ~(isnumeric(phi_p) && phi_p >= 0 && phi_p < 1)
    error('phi_p must be a number in [0,1)')
end
if ~(isnumeric(Theta) && Theta >= 0)
    error('Theta must be a nonnegative number')
end
if ~(isnumeric(L) && L >= 0)
    error('L must be a nonnegative number')
end

% define simplifying variables
R_m = (((3+phi_m-4*phi_m^2)/3)*phi_m^(1+2*Theta))^(-1/2);
R_p = (((3+phi_p-4*phi_p^2)/3)*phi_p^(1+2*Theta))^(-1/2);
E_m = exp(-R_m*L);
E_p = exp(-R_p*L);
F_m = phi_m^(2+2*Theta)*(1-phi_m);
F_p = phi_p^(2+2*Theta)*(1-phi_p);
S_m = ((phi_m-4*phi_m^2)/(3+phi_m-4*phi_m^2));
S_p = ((phi_p-4*phi_p^2)/(3+phi_p-4*phi_p^2));

% check that phi_m and phi_p are large enough
if ~isinf(exp(R_m*L)) && ~isinf(exp(R_p*L))
    
    % solve linear system for unknown constants
    A = [E_p                  1                   0                  0        0  0 ; ...
        0                     0                   1                 E_m       0  0 ; ...
        1/R_p              -E_p/R_p               0                  0        1  0 ; ...
        0                     0                E_m/R_m            -1/R_m      0  1 ; ...
        F_p                 F_p*E_p           -F_m*E_m             -F_m       0  0 ; ...
        -phi_m*F_p*R_p  phi_m*F_p*R_p*E_p  phi_p*F_m*R_m*E_m  -phi_p*F_m*R_m  0  0 ];
    
    b = [1 1 0 0 F_p-F_m 0]';
    
    c = A\b; C1_p = c(1); C2_p = c(2); C1_m = c(3); C2_m = c(4); CInt_p = c(5);  CInt_m = c(6);
    
    % if phi_m is very close to zero
elseif isinf(exp(R_m*L)) && ~isinf(exp(R_p*L))
    
    % solve simplified system for unknown constants
    C1_p = 1/(1+E_p);
    C2_p = 1/(1+E_p);
    CInt_p = -(1-E_p)/((1+E_p)*R_p);
    C1_m = 0; C2_m = 0; CInt_m = 0;
    
    % if phi_m is very close to zero
elseif ~isinf(exp(R_m*L)) && isinf(exp(R_p*L))
    
    % solve simplified system for unknown constants
    C1_m = 1/(1+E_m);
    C2_m = 1/(1+E_m);
    CInt_m = -(1-E_m)/((1+E_m)*R_m);
    C1_p = 0; C2_p = 0; CInt_p = 0;
end

% avoid multiplying inf by 0
exp_m = @(x) min(exp(R_m*x),1);
exp_p = @(x) min(exp(-R_p*x),1);

%   1. Darcy or scaled velocity
if FORMULATION == 1 || FORMULATION == 2
    field(1).uExact{1} = @(x,y)                  -F_m*(1 - C1_m*exp(-R_m*(x+L)) - C2_m*exp_m(x))*(x<0) - ...
        F_p*(1 - C1_p*exp_p(x) - C2_p*exp(R_p*(x-L)))*(x>=0);
else
    field(1).uExact{1} = @(x,y) -phi_m^(-1-Theta)*F_m*(1 - C1_m*exp(-R_m*(x+L)) - C2_m*exp_m(x))*(x<0) - ...
        phi_p^(-1-Theta)*F_p*(1 - C1_p*exp_p(x) - C2_p*exp(R_p*(x-L)))*(x>=0);
end
if FORMULATION == 1 || FORMULATION == 2
    field(1).gradUExact{1,1} = @(x,y)                  -F_m*(C1_m*R_m*exp(-R_m*(x+L)) - C2_m*R_m*exp_m(x))*(x<0) - ...
        F_p*(C1_p*R_p*exp_p(x) - C2_p*R_p*exp(R_p*(x-L)))*(x>=0);
else
    field(1).gradUExact{1,1} = @(x,y) -phi_m^(-1-Theta)*F_m*(C1_m*R_m*exp(-R_m*(x+L)) - C2_m*R_m*exp_m(x))*(x<0) - ...
        phi_p^(-1-Theta)*F_p*(C1_p*R_p*exp_p(x) - C2_p*R_p*exp(R_p*(x-L)))*(x>=0);
end

% 2. fluid or scaled fluid pressure potential
if FORMULATION == 4
    field(2).uExact{1} = @(x,y) sqrt(phi_m)*(1-phi_m)*(x + C1_m*(1/R_m)*exp(-R_m*(x+L)) - C2_m*(1/R_m)*exp_m(x) + CInt_m)*(x<0) + ...
        sqrt(phi_p)*(1-phi_p)*(x + C1_p*(1/R_p)*exp_p(x) - C2_p*(1/R_p)*exp(R_p*(x-L)) + CInt_p)*(x>=0);
else
    field(2).uExact{1} = @(x,y)             (1-phi_m)*(x + C1_m*(1/R_m)*exp(-R_m*(x+L)) - C2_m*(1/R_m)*exp_m(x) + CInt_m)*(x<0) + ...
        (1-phi_p)*(x + C1_p*(1/R_p)*exp_p(x) - C2_p*(1/R_p)*exp(R_p*(x-L)) + CInt_p)*(x>=0);
end

% 3. matrix velocity
field(3).uExact{1} = @(x,y) F_m*(1 - C1_m*exp(-R_m*(x+L)) - C2_m*exp_m(x))*(x<0) + ...
    F_p*(1 - C1_p*exp_p(x) - C2_p*exp(R_p*(x-L)))*(x>=0);
field(3).gradUExact{1,1} = @(x,y) F_m*(C1_m*R_m*exp(-R_m*(x+L)) - C2_m*R_m*exp_m(x))*(x<0) + ...
    F_p*(C1_p*R_p*exp_p(x) - C2_p*R_p*exp(R_p*(x-L)))*(x>=0);

% 4. matrix prssure potential
field(4).uExact{1} = @(x,y) (1-phi_m)*(x + C1_m*(S_m/R_m)*exp(-R_m*(x+L)) - C2_m*(S_m/R_m)*exp_m(x) + CInt_m)*(x<0) + ...
    (1-phi_p)*(x + C1_p*(S_p/R_p)*exp_p(x) - C2_p*(S_p/R_p)*exp(R_p*(x-L)) + CInt_p)*(x>=0);

% 5.auxilliary velocity (if exists)
if FORMULATION == 2
    field(5).uExact{1} = @(x,y) -(1-phi_m)*(1 - C1_m*exp(-R_m*(x+L)) - C2_m*exp(R_m*x))*(x<0) - ...
        (1-phi_p)*(1 - C1_p*exp(-R_p*x) - C2_p*exp(R_p*(x-L)))*(x>=0);
    field(5).gradUExact{1,1} = @(x,y) (1-phi_m)*(C1_m*R_m*exp(-R_m*(x+L)) - C2_m*R_m*exp(R_m*x))*(x<0) + ...
        (1-phi_p)*(C1_p*R_p*exp(-R_p*x) - C2_p*R_p*exp(R_p*(x-L)))*(x>=0);
end

field(1).gradUExact{1,2} = @(x,y) 0;
field(1).uExact{2} = @(x,y) 0;
field(1).gradUExact{2,1} = @(x,y) 0;
field(1).gradUExact{2,2} = @(x,y) 0;

field(2).gradUExact{1,1} = @(x,y) 0;
field(2).gradUExact{1,2} = @(x,y) 0;

field(3).gradUExact{1,2} = @(x,y) 0;
field(3).uExact{2} = @(x,y) 0;
field(3).gradUExact{2,1} = @(x,y) 0;
field(3).gradUExact{2,2} = @(x,y) 0;

field(4).gradUExact{1,1} = @(x,y) 0;
field(4).gradUExact{1,2} = @(x,y) 0;

if FORMULATION == 2
    field(5).gradUExact{1,2} = @(x,y) 0;
    field(5).uExact{2} = @(x,y) 0;
    field(5).gradUExact{2,1} = @(x,y) 0;
    field(5).gradUExact{2,2} = @(x,y) 0;
end