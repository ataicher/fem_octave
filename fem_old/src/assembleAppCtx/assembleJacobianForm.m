%% JACOBAIN FORM-----------------------------------------------------------
% given variational form,
%
% F(u) = (v , f0(u,gradU)) + (gradV , f1(u,gradU))
%
% let,
%
% u = sum u_k phi_k,
% gradU = sum u_k gradPhi_k
%
% let phi be a test function in field fI and psi a test function in field
% fJ. Note that with a vector field, each term must be decomposed into
% scalar products and then summed using Einstien summation notation,
% summing over the components i, j and derivative directions k,l
%
% dF/du = (df0/du(u,gradU,x,y)_i,j * phi_i, psi_j ) + ...
%         (df0/d(gradU)(u,gradU)_i,kj * gradPhi_ik, psi_j ) + ...
%         (df1/du(u,gradU)_i,jl * phi_i, gradPsi_jl ) + ...
%         (df1/d(gradU)(u,gradU)_ik,jl * gradPhi_ik, gradPsi_jl )
%
% INPUT:
%   df00:   df00{fI,fJ}{i,j} = df0/du(u,gradU)_i,j
%   df01:   df01{fI,fJ}{i,j,l} = df0/d(gradU)(u,gradU)_i,jl
%   df10:   df10{fI,fJ}{i,k,j} = df1/du(u,gradU)_ik,j
%   df11:   df11{fI,fJ}{i,k,j,l} = df1/d(gradU)(u,gradU)_ik,jl
%   h:      finite difference parameter.  if h = 0, then use user defined
%           Jacobian form. if h > 1 then use a finite difference scheme to
%           approximate the Jacobian using a finite difference parameter h
%   appCtx: struct containg problem information
%
% OUTPUT: jacFrom struct containing
%   df00:   same as input
%   df01:   same as input
%   df10:   same as input
%   df11:   same as input
%
function appCtx = assembleJacobianForm(appCtx)

dim = appCtx.dim;
numFields = appCtx.numFields;

for fI=1:numFields
    JacForm = appCtx.field(fI).JacForm;
    numCompI = appCtx.field(fI).numComp;
    
    if length(JacForm.field) < numFields
        JacForm.field(numFields).v = [];
    end
    
    for fJ=1:numFields
        numCompJ = appCtx.field(fJ).numComp;
        
        if ~isfield(JacForm.field(fJ), 'v') || length(JacForm.field(fJ).v) < numCompI
            JacForm.field(fJ).v{numCompI} = [];
        end
        if ~isfield(JacForm.field(fJ), 'gradV') || size(JacForm.field(fJ).gradV,1) < numCompI ||  size(JacForm.field(fJ).gradV,2) < dim
            JacForm.field(fJ).gradV{numCompI,dim} = [];
        end
        
        for compI = 1:numCompI
            if ~isfield(JacForm.field(fJ).v{compI}, 'w') || length(JacForm.field(fJ).v{compI}.w) < numCompJ
                JacForm.field(fJ).v{compI}.w{numCompJ} = [];
            end
            if ~isfield(JacForm.field(fJ).v{compI}, 'gradW') || size(JacForm.field(fJ).v{compI}.gradW,1) < numCompJ ||  size(JacForm.field(fJ).v{compI}.gradW,2) < dim
                JacForm.field(fJ).v{compI}.gradW{numCompJ,dim} = [];
            end
            
            for d=1:dim
                if ~isfield(JacForm.field(fJ).gradV{compI,d}, 'w') || length(JacForm.field(fJ).gradV{compI,d}.w) < numCompJ
                    JacForm.field(fJ).gradV{compI,d}.w{numCompJ} = [];
                end
                if ~isfield(JacForm.field(fJ).gradV{compI,d}, 'gradW') || size(JacForm.field(fJ).gradV{compI,d}.gradW,1) < numCompJ ||  size(JacForm.field(fJ).gradV{compI,d}.gradW,2) < dim
                    JacForm.field(fJ).gradV{compI,d}.gradW{numCompJ,dim} = [];
                end
            end
        end 
    end
    
    appCtx.field(fI).JacForm = JacForm;
end

% unset terms set to zero
for fI=1:numFields
    JacForm = appCtx.field(fI).JacForm;
    numCompI = appCtx.field(fI).numComp;
    
    for fJ=1:numFields
        numCompJ = appCtx.field(fJ).numComp;
        
        for compI = 1:numCompI
            
            for compJ = 1:numCompJ
                
                if isempty(JacForm.field(fJ).v{compI}.w{compJ})
                    JacForm.field(fJ).v{compI}.w{compJ} = @(u,gradU,x,y) 0;
                else
                    JacForm.field(fJ).v{compI}.w{compJ} = expandFun(JacForm.field(fJ).v{compI}.w{compJ});
                end
                
                for d=1:dim
                    if isempty(JacForm.field(fJ).v{compI}.gradW{compJ,d})
                        JacForm.field(fJ).v{compI}.gradW{compJ,d} = @(u,gradU,x,y) 0;
                    else
                        JacForm.field(fJ).v{compI}.gradW{compJ,d} = expandFun(JacForm.field(fJ).v{compI}.gradW{compJ,d});
                    end
                    if isempty(JacForm.field(fJ).gradV{compI,d}.w{compJ})
                        JacForm.field(fJ).gradV{compI,d}.w{compJ} = @(u,gradU,x,y) 0;
                    else
                        JacForm.field(fJ).gradV{compI,d}.w{compJ} = expandFun(JacForm.field(fJ).gradV{compI,d}.w{compJ});
                    end
                    
                    for d2=1:dim
                        if isempty(JacForm.field(fJ).gradV{compI,d}.gradW{compJ,d2})
                            JacForm.field(fJ).gradV{compI,d}.gradW{compJ,d2} = @(u,gradU,x,y) 0;
                        else
                            JacForm.field(fJ).gradV{compI,d}.gradW{compJ,d2} = expandFun(JacForm.field(fJ).gradV{compI,d}.gradW{compJ,d2});
                        end
                    end
                end
            end
        end
    end
    
    appCtx.field(fI).JacForm = JacForm;
end

% convert to strings
for fI=1:numFields
    JacForm = appCtx.field(fI).JacForm;
    numCompI = appCtx.field(fI).numComp;
    
    for fJ=1:numFields
        numCompJ = appCtx.field(fJ).numComp;
        
        for compI = 1:numCompI
            
            for compJ = 1:numCompJ
                str = func2str(JacForm.field(fJ).v{compI}.w{compJ});
                JacFormStr.field(fJ).v{compI}.w{compJ} = str(15:end);
                
                for d=1:dim
                    str = func2str(JacForm.field(fJ).v{compI}.gradW{compJ,d});
                    JacFormStr.field(fJ).v{compI}.gradW{compJ,d} = str(15:end);
                    str = func2str(JacForm.field(fJ).gradV{compI,d}.w{compJ});
                    JacFormStr.field(fJ).gradV{compI,d}.w{compJ} = str(15:end);
                    
                    for d2=1:dim
                        str = func2str(JacForm.field(fJ).gradV{compI,d}.gradW{compJ,d2});
                        JacFormStr.field(fJ).gradV{compI,d}.gradW{compJ,d2} = str(15:end);
                    end
                end
            end
        end
    end
    
    appCtx.field(fI).JacFormStr =  JacFormStr;
end

% convert to vector valued functions
appCtx.field = rmfield(appCtx.field, 'JacForm');
for fI=1:numFields
    JacFormStr = appCtx.field(fI).JacFormStr;
    numCompI = appCtx.field(fI).numComp;
    
    for fJ=1:numFields
        numCompJ = appCtx.field(fJ).numComp;
        
        str = '@(u,gradU,x,y)[ ';
        for compI = 1:numCompI
            for compJ = 1:numCompJ
                str = [ str , JacFormStr.field(fJ).v{compI}.w{compJ} , '; '];
            end
        end
        str = [str(1:end-2), ']'];
        appCtx.field(fI).JacForm.field(fJ).VW = str2func(str);
        
        str = '@(u,gradU,x,y)[ ';
        for compI = 1:numCompI
            for compJ = 1:numCompJ
                for d=1:dim
                    str = [ str , JacFormStr.field(fJ).v{compI}.gradW{compJ,d}, '; '];
                end
            end
        end
        str = [str(1:end-2), ']'];
        appCtx.field(fI).JacForm.field(fJ).VGradW = str2func(str);
        
        str = '@(u,gradU,x,y)[ ';
        for compI = 1:numCompI
            for d=1:dim
                for compJ = 1:numCompJ
                    str = [ str , JacFormStr.field(fJ).gradV{compI,d}.w{compJ}, '; '];
                end
            end
        end
        str = [str(1:end-2), ']'];
        appCtx.field(fI).JacForm.field(fJ).gradVW = str2func(str);
        
        str = '@(u,gradU,x,y)[ ';
        for compI = 1:numCompI
            for d=1:dim
                for compJ = 1:numCompJ
                    for d2=1:dim
                        str = [ str , JacFormStr.field(fJ).gradV{compI,d}.gradW{compJ,d2}, '; '];
                    end
                end
                
            end
        end
        str = [str(1:end-2), ']'];
        appCtx.field(fI).JacForm.field(fJ).gradVGradW = str2func(str);
    end
end

appCtx.field = rmfield(appCtx.field, 'JacFormStr');