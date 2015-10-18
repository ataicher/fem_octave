% Compute residual. The variational residual statement gives for each field
% an equation:
%   <v,f(u)> = int(v * f0(u,gradU) + gradV : f1(u,gradU)) = 0
% Can be writen as:
%   sum_e E_e' * [ (B_e)'*W_e * Lambda(f_0(u_e,gradU_e))
%       + sum_d (D_e)' * W_e * lambda(f_1(u_e,gradU_e)) ]= 0
%   where,
%       u_e = B_e * E_e * u
%       (gradU_e)_i = (D_e)_i * E_e * u
%
% INPUT:
%   u - global degrees of freedom
%   appCtx - problem parameters
% OUTPUT:
%   r - residual
%
function Connectivity = computeConnectivityMatrix(appCtx)
LG = appCtx.LG;
variables = appCtx.variables;
mesh = appCtx.mesh;
FE = appCtx.FE;
numFields = variables.numFields;
fNames = variables.fNames;
numCells = mesh.numCells;


% initialize Connectivity matrix
Connectivity = zeros(globalSize(appCtx),globalSize(appCtx));
for c=1:numCells
    
    for fI=1:numFields
        fNameI = fNames{fI};
        numCompI = variables.numComp.(fNameI);
        numBasisFuncsI = FE.numBasisFuncs.(fNameI);
        
        for fJ=1:numFields
            fNameJ = fNames{fJ};
            numCompJ = variables.numComp.(fNameJ);
            numBasisFuncsJ = FE.numBasisFuncs.(fNameJ);
            
            for compI=1:numCompI
                for bI=1:numBasisFuncsI
                    for compJ=1:numCompJ
                        for bJ=1:numBasisFuncsJ
                            Connectivity(LG.cell.(fNameI)(bI,c),LG.cell.(fNameJ)(bJ,c)) = 1;
                        end
                    end
                end
            end
        end
    end
end