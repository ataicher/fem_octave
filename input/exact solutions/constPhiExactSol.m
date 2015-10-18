if ~(exist('phi_0', 'var') && exist('Theta', 'var') && exist('L', 'var'))
    error('phi_0, Theta, or L does not exist')
end
if ~(isnumeric(phi_0) && phi_0 >= 0 && phi_0 < 1)
    error('phi_0 must be a number in [0,1)')
end
if ~(isnumeric(Theta) && Theta >= 0)
    error('Theta must be a nonnegative number')
end
if ~(isnumeric(L) && L >= 0)
    error('L must be a nonnegative number')
end

R = (((3+phi_0-4*phi_0^2)/3)*phi_0^(1+2*Theta))^(-1/2);
F = phi_0^(2+2*Theta)*(1-phi_0);
S = (phi_0-4*phi_0^2)/(3+phi_0-4*phi_0^2);

if FORMULATION == 3 || FORMULATION == 4
    field(1).uExact{1} = @(x,y) -phi_0^(-1-Theta)*F*(1 - cosh(R*x)/cosh(R*L));
    field(1).gradUExact{1,1} = @(x,y) phi_0^(-1-Theta)*F*R*sinh(R*x)/cosh(R*L);
else
    field(1).uExact{1} = @(x,y)                  -F*(1 - cosh(R*x)/cosh(R*L));
    field(1).gradUExact{1,1} = @(x,y)                  F*R*sinh(R*x)/cosh(R*L);
end

if FORMULATION == 4
    field(2).uExact{1} = @(x,y) sqrt(phi_0)*(1-phi_0)*(x - (1/R)*sinh(R*x)/cosh(R*L));
else
    field(2).uExact{1} = @(x,y)             (1-phi_0)*(x - (1/R)*sinh(R*x)/cosh(R*L));
end

field(3).uExact{1} = @(x,y) F*(1 - cosh(R*x)/cosh(R*L));
field(3).gradUExact{1,1} = @(x,y) -F*R*sinh(R*x)/cosh(R*L);

field(4).uExact{1} = @(x,y) (1-phi_0)*(x - (S/R)*sinh(R*x)/cosh(R*L));

if FORMULATION == 2
    field(5).uExact{1} = @(x,y) -phi_0^(-2-2*Theta)*F*(1 - cosh(R*x)/cosh(R*L));
    field(5).gradUExact{1,1} = @(x,y) phi_0^(-2-2*Theta)*F*R*sinh(R*x)/cosh(R*L);
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