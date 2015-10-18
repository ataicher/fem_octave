% generate exact solution for constant porosity for 1D
% compacting column with no flow BC on [-L,L]
%       phi = phi_0

R = (((3-phi_0+4*phi_0^2)/3)*phi_0^(1+2*Theta))^(-1/2);
C = 1/(1+exp(-2*R*L));

if SETDEGENERATE || NAIVE == 2
    field(1).uExact{1} = @(x,y) -(1-phi_0)*phi_0^(1+Theta)*(1 - C*exp(-R*(x+L)) - C*exp(R*(x-L)));
    field(1).gradUExact{1,1} = @(x,y) -(1-phi_0)*phi_0^(1+Theta)*C*R*(exp(-R*(x+L)) - exp(R*(x-L)));
else
    field(1).uExact{1} = @(x,y) -(1-phi_0)*phi_0^(2+2*Theta)*(1 - C*exp(-R*(x+L)) - C*exp(R*(x-L)));
    field(1).gradUExact{1,1} = @(x,y) -(1-phi_0)*phi_0^(2+2*Theta)*C*R*(exp(-R*(x+L)) - exp(R*(x-L)));
end
field(1).gradUExact{1,2} = @(x,y) 0;
field(1).uExact{2} = @(x,y) 0;
field(1).gradUExact{2,1} = @(x,y) 0;
field(1).gradUExact{2,2} = @(x,y) 0;

if SETDEGENERATE
    field(2).uExact{1} = @(x,y) sqrt(phi_0)*(1-phi_0)*(x + (C/R)*(exp(-R*(x+L)) - exp(R*(x-L))));
    field(2).gradUExact{1,1} = @(x,y) 0;
    field(2).gradUExact{1,2} = @(x,y) 0;
else
    field(2).uExact{1} = @(x,y) (1-phi_0)*(x + (C/R)*(exp(-R*(x+L)) - exp(R*(x-L))));
    field(2).gradUExact{1,1} = @(x,y) 0;
    field(2).gradUExact{1,2} = @(x,y) 0;
end

field(3).uExact{1} = @(x,y) (1-phi_0)*phi_0^(2+2*Theta)*(1 - C*exp(-R*(x+L)) - C*exp(R*(x-L)));
field(3).gradUExact{1,1} = @(x,y) (1-phi_0)*phi_0^(2+2*Theta)*C*R*(exp(-R*(x+L)) - exp(R*(x-L)));
field(3).gradUExact{1,2} = @(x,y) 0;
field(3).uExact{2} = @(x,y) 0;
field(3).gradUExact{2,1} = @(x,y) 0;
field(3).gradUExact{2,2} = @(x,y) 0;

field(4).uExact{1} = @(x,y) (1-phi_0)*(x - ((1-4*phi_0)/(3-phi_0+4*phi_0^2))*phi_0*(C/R)*(exp(-R*(x+L)) - exp(R*(x-L))));
field(4).gradUExact{1,1} = @(x,y) 0;
field(4).gradUExact{1,2} = @(x,y) 0;

if NAIVE == 4 && ~SETDEGENERATE
    field(5).uExact{1} = @(x,y) -(1-phi_0)*(1 - C*exp(-R*(x+L)) - C*exp(R*(x-L)));
    field(5).gradUExact{1,1} = @(x,y) -(1-phi_0)*C*R*(exp(-R*(x+L)) - exp(R*(x-L)));
    field(5).gradUExact{1,2} = @(x,y) 0;
    field(5).uExact{2} = @(x,y) 0;
    field(5).gradUExact{2,1} = @(x,y) 0;
    field(5).gradUExact{2,2} = @(x,y) 0;
end