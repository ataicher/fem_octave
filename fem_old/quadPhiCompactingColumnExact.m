% generate exact solution approximation for quadratically increasing porosity for 1D
% compacting column with no flow BC on [-L,L]
%             { phi_m               x <= 0
%       phi = { phi_m + C*x^2   0 < x <= ((phi_p-phi_m)/C)^(1/2)
%             { phi_p               x > ((phi_p-phi_m)/C)^(1/2)
% only holds for Theta = 0

r1 = (-3 + sqrt(13))/2;
r2 = (-3 - sqrt(13))/2;

if SETDEGENERATE || NAIVE == 2
    field(1).uExact{1} = @(x,y) (x^2/((1-r1)*(1-r2)))*(x^4 - L^(1-r1)*x^(r1+3))*(x>0);
    field(1).gradUExact{1,1} = @(x,y) (1/((1-r1)*(1-r2)))*(4*x^3 - (r1+3)*L^(1-r1)*x^(r1+2))*(x>0);
else
    field(1).uExact{1} = @(x,y) (1/((1-r1)*(1-r2)))*(x^4 - L^(1-r1)*x^(r1+3))*(x>0);
    field(1).gradUExact{1,1} = @(x,y) (x^2/((1-r1)*(1-r2)))*(4*x^3 - (r1+3)*L^(1-r1)*x^(r1+2))*(x>0);
end
field(1).gradUExact{1,2} = @(x,y) 0;
field(1).uExact{2} = @(x,y) 0;
field(1).gradUExact{2,1} = @(x,y) 0;
field(1).gradUExact{2,2} = @(x,y) 0;

if SETDEGENERATE
    field(2).uExact{1} = @(x,y) (x/(r1*(1-r1)*(1-r2)))*(-r1*x + L^(1-r1)*x^r1)*(x>0);
    field(2).gradUExact{1,1} = @(x,y) 0;
    field(2).gradUExact{1,2} = @(x,y) 0;
else
    field(2).uExact{1} = @(x,y) (1/(r1*(1-r1)*(1-r2)))*(-r1*x + L^(1-r1)*x^r1)*(x>0);
    field(2).gradUExact{1,1} = @(x,y) 0;
    field(2).gradUExact{1,2} = @(x,y) 0;
end
field(3).uExact{1} = @(x,y) -(1/((1-r1)*(1-r2)))*(x^4 - L^(1-r1)*x^(r1+3))*(x>0);
field(3).gradUExact{1,1} = @(x,y) -(1/((1-r1)*(1-r2)))*(4*x^3 - (r1+3)*L^(1-r1)*x^(r1+2))*(x>0);
field(3).gradUExact{1,2} = @(x,y) 0;
field(3).uExact{2} = @(x,y) 0;
field(3).gradUExact{2,1} = @(x,y) 0;
field(3).gradUExact{2,2} = @(x,y) 0;

field(4).uExact{1} = @(x,y) x;
field(4).gradUExact{1,1} = @(x,y) 0;
field(4).gradUExact{1,2} = @(x,y) 0;

if NAIVE == 4 && ~SETDEGENERATE
    field(1).uExact{1} = @(x,y) (x^-2/((1-r1)*(1-r2)))*(x^4 - L^(1-r1)*x^(r1+3))*(x>0);
    field(1).gradUExact{1,1} = @(x,y) (x^-2/((1-r1)*(1-r2)))*(4*x^3 - (r1+3)*L^(1-r1)*x^(r1+2))*(x>0);
    field(5).gradUExact{1,1} = @(x,y) 0;
    field(5).gradUExact{1,2} = @(x,y) 0;
    field(5).uExact{2} = @(x,y) 0;
    field(5).gradUExact{2,1} = @(x,y) 0;
    field(5).gradUExact{2,2} = @(x,y) 0;
end
