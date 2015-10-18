if ~(exist('phi', 'var') && exist('Theta', 'var') && exist('L', 'var'))
    error('phi, Theta, or L does not exist')
end
if ~isa(phi, 'function_handle')
    error('phi must be a function handle')
end
if ~(isnumeric(Theta) && Theta >= 0)
    error('Theta must be a nonnegative number')
end
if ~(isnumeric(L) && L >= 0)
    error('L must be a nonnegative number')
end
r1 = (-3 + sqrt(13))/2;
r2 = (-3 - sqrt(13))/2;

if FORMULATION == 3 || FORMULATION == 4
    field(1).uExact{1} = @(x,y) (1/((1-r1)*(1-r2)))*(x^2 - L^(1-r1)*x^(r1+1))*(x>0);
    field(1).gradUExact{1,1} = @(x,y) (1/((1-r1)*(1-r2)))*(2*x - (r1+1)*L^(1-r1)*x^r1)*(x>0);
else
    field(1).uExact{1} = @(x,y) (1/((1-r1)*(1-r2)))*(x^4 - L^(1-r1)*x^(r1+3))*(x>0);
    field(1).gradUExact{1,1} = @(x,y) (1/((1-r1)*(1-r2)))*(4*x^3 - (r1+3)*L^(1-r1)*x^(r1+2))*(x>0);
end

if FORMULATION == 4
    field(2).uExact{1} = @(x,y) (1/(r1*(1-r1)*(1-r2)))*(-r1*x^2 + L^(1-r1)*x^(r1+1))*(x>0);
else
    field(2).uExact{1} = @(x,y) (1/(r1*(1-r2))) + (1/(r1*(1-r1)*(1-r2)))*(-r1*x + L^(1-r1)*x^r1)*(x>0);
end

field(3).uExact{1} = @(x,y) (1/((1-r1)*(1-r2)))*(-x^4 + L^(1-r1)*x^(r1+3))*(x>0);
field(3).gradUExact{1,1} = @(x,y) (1/((1-r1)*(1-r2)))*(-4*x^3 + (r1+3)*L^(1-r1)*x^(r1+2))*(x>0);

if FORMULATION == 4
    field(4).uExact{1} = @(x,y) x;
else
    field(4).uExact{1} = @(x,y) (1/(r1*(1-r2))) + x;
end

if FORMULATION == 2
    field(5).uExact{1} = @(x,y) (1/((1-r1)*(1-r2)))*(1 - L^(1-r1)*x^(r1-1))*(x>0);
    field(5).gradUExact{1,1} = @(x,y) (1/(1-r2))*(L^(1-r1)*x^(r1-2))*(x>0);
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