%% INITIALIZE EXACT SOLUTION
%   initialize function handles for the exact solution and its first and
%   second derivatives
% INPUT:
%   variables: variable information
%   fun: struct contatining strings with the exact solution
% OUTPUT:
%   for each field f and component c,
%       fun: exact solution
%       fun_X: exact solution x-derivative
%       fun_Y: exact solution y-derivative
%       fun_XX: exact solution 2nd derivative
%       fun_XY: exact solution 2nd derivative
%       fun_YX: exact solution 2nd derivative
%       fun_YY: exact solution 2nd derivative
%
function field = symbolicDerivatives(field, order)

numFields = length(field);
syms x y;

if length(order) == 1
    order = order*ones(1,numFields);
end

for f = 1:numFields;
    numComp = field(f).numComp;
    
    for comp = 1:numComp
        
        if order(f) == 1
            funTmp = matlabFunction(field(f).uExact{comp},'vars',[x,y]);
            field(f).gradUExact{comp,1} = matlabFunction(diff(funTmp(x,y),x),'vars',[x,y]);
            field(f).gradUExact{comp,2} = matlabFunction(diff(funTmp(x,y),y),'vars',[x,y]);
        elseif order(f) == 2
            funTmp = matlabFunction(field(f).uExact{comp},'vars',[x,y]);
            field(f).gradUExact{comp,1} = matlabFunction(diff(funTmp(x,y),x),'vars',[x,y]);
            field(f).gradUExact{comp,2} = matlabFunction(diff(funTmp(x,y),y),'vars',[x,y]);
            field(f).grad2UExact{comp,1,1} = matlabFunction(diff(funTmp(x,y),x,2),'vars',[x,y]);
            field(f).grad2UExact{comp,1,2} = matlabFunction(diff(funTmp(x,y),x,y),'vars',[x,y]);
            field(f).grad2UExact{comp,2,1} = matlabFunction(diff(funTmp(x,y),y,x),'vars',[x,y]);
            field(f).grad2UExact{comp,2,2} = matlabFunction(diff(funTmp(x,y),y,2),'vars',[x,y]);
        else
            error('myApp:argChk', 'only differentiates up to second derivatives')
        end
    end
end