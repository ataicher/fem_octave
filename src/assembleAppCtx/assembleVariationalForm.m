%% ASSEMBLE VARIATIONAL FORM
% assemble variational form information.
% rewriting the variational form, for all v in V
%
% F(u) = (f0(u,gradU), v) + (f1(u,gradU), gradV)
%
% note that with a vector field, each term must be decomposed into scalar
% products and then summed using Einstien summation notation, summing over 
% the i components and j derivative directions
%
% F(u) = (f0(u,gradU)_i, v_i) + (f1(u,gradU)_ij, gradV_ij)
%
% INPUT:
%   f0:     function handle representing part of variational form 
%           multiplied by test function for field f
%           f0{f}{comp}(u,gradU,x,y) - function value at domain point (x,y)
%           for variational form of component comp of field f
%   f1:     function handle representing part of variational form
%           multiplied by test function derivative of field f
%           f1{f}{comp,d}(u,gradU,x,y) - function value at domain point 
%           (x,y) for variational form of d-derivative of component comp of
%           field f
%   appCtx: struct containg problem information
%
% OUTPUT: varForm struct containing
%   f0: same as input
%   f1: same as input
%
function appCtx = assembleVariationalForm(appCtx)

dim = appCtx.dim;
numFields = appCtx.numFields;

% resize if necessary
for f=1:numFields
    numComp = appCtx.field(f).numComp;
    
    if ~isfield(appCtx.field(f).varForm, 'v') || length(appCtx.field(f).varForm.v) < numComp
        appCtx.field(f).varForm.v{numComp} = [];
    end
    if ~isfield(appCtx.field(f).varForm, 'gradV') || size(appCtx.field(f).varForm.gradV,1) < numComp ||  size(appCtx.field(f).varForm.gradV,2) < dim
        appCtx.field(f).varForm.gradV{numComp,dim} = [];
    end
end

% unset terms set to zero and expand functions
for f=1:numFields
    numComp = appCtx.field(f).numComp;   
    varForm = appCtx.field(f).varForm;
    
    for comp = 1:numComp      
        if isempty(varForm.v{comp})
            varForm.v{comp} = @(u,gradU,x,y) 0;
        else
            %varForm.v{comp} = expandFun(varForm.v{comp});
        end       
        for d=1:dim            
            if isempty(varForm.gradV{comp,d})
                varForm.gradV{comp,d} = @(u,gradU,x,y) 0;
            else
                %varForm.gradV{comp,d} = expandFun(varForm.gradV{comp,d});
            end           
        end
    end
    
    appCtx.field(f).varForm = varForm;
end
