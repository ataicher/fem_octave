% calculate the L2 error from the degrees of freedom using values at
% quadrature points
% INPUT:
%   u       - degrees of freedoms
%   appCtx  - problem specific information
% OUTPUT:
%   error   - error(field) is the H1 error of the field
function [norms, varargout] = calcNormAndRelError(u, appCtx)

dim = appCtx.dim;
numCells = appCtx.mesh.numCells;
numQuadPoints = appCtx.quad.numQuadPoints;
quadWeights = appCtx.quad.quadWeights;
numFields = appCtx.numFields;
EXISTEXACTSOL = appCtx.EXISTEXACTSOL;
EXISTEXACTSOLGRAD = appCtx.EXISTEXACTSOLGRAD;

if EXISTEXACTSOL
    uExactDOF = projectExactSolution(appCtx);
end

% initialize norm and error structs
for f = 1:numFields
    norms(f).L2 = 0;
    norms(f).L2Exact = 0;
    norms(f).L2ExactInt = 0;
    if EXISTEXACTSOL
        norms(f).L2Exact = 0;
        norms(f).L2ExactInt = 0;
        error(f).L2 = 0;
        error(f).L2Int = 0;
    end
    if appCtx.field(f).CALCHDIVNORM
        norms(f).HDiv = 0;
        if EXISTEXACTSOLGRAD
            norms(f).HDivExact = 0;
            norms(f).HDivExactInt = 0;
            error(f).HDiv = 0;
            error(f).HDivInt = 0;
        end
    end
    if appCtx.field(f).CALCH1NORM
        norms(f).H1 = 0;
        if EXISTEXACTSOLGRAD
            norms(f).H1Exact = 0;
            norms(f).H1ExactInt = 0;
            error(f).H1 = 0;
            error(f).H1Int = 0;
        end
    end
end

for c=1:numCells
    detJ = appCtx.cellGeometry.detJ{c};
    realBasis = projectBasis(c,appCtx);
    realBasisDer = projectBasisDer(c,appCtx);
    realQuadPoints = projectQuadPoints(c,appCtx);
    
    for f=1:numFields
        numComp = appCtx.field(f).numComp;
        cellDOF = appCtx.field(f).cellDOF(:,c);
        rBasis = realBasis.field{f};
        rBasisDer = realBasisDer.field{f};
        
        if EXISTEXACTSOL
            uExact = appCtx.field(f).uExact;
            uExactLocal = uExactDOF(cellDOF);
            if EXISTEXACTSOLGRAD
                gradUExact = appCtx.field(f).gradUExact;
            end
        end
        uLocal = u(cellDOF);
        
        for q=1:numQuadPoints
            for comp=1:numComp
                uVal = uLocal'*rBasis(:,comp,q);
                gradUVal = uLocal'*rBasisDer(:,:,comp,q);
                
                % L2 norm
                norms(f).L2 = norms(f).L2 + uVal^2*detJ*quadWeights(q);
                
                if EXISTEXACTSOL
                    uExactVal = uExact{comp}(realQuadPoints(1,q),realQuadPoints(2,q));
                    uExactIntVal = uExactLocal'*rBasis(:,comp,q);
                    
                    % L2 norm exact
                    norms(f).L2Exact = norms(f).L2Exact + uExactIntVal^2*detJ*quadWeights(q);
                    
                    % L2 interpolation norm exact
                    norms(f).L2ExactInt = norms(f).L2ExactInt + uExactIntVal^2*detJ*quadWeights(q);
                    
                    %L2 error
                    error(f).L2 = error(f).L2 + (uVal-uExactVal)^2*detJ*quadWeights(q);
                    
                    % L2 interpolation error
                    error(f).L2Int = error(f).L2Int + (uVal-uExactIntVal)^2*detJ*quadWeights(q);                    
                end
                
                if EXISTEXACTSOLGRAD
                    gradUExactVal = [gradUExact{comp,1}(realQuadPoints(1,q),realQuadPoints(2,q)); ...
                        gradUExact{comp,2}(realQuadPoints(1,q),realQuadPoints(2,q))];
                    gradUExactIntVal = uExactLocal'*rBasisDer(:,:,comp,q);
                end
                
                if appCtx.field(f).CALCHDIVNORM
                    
                    % H(div) norm
                    norms(f).HDiv = norms(f).HDiv + gradUVal(comp)^2*detJ*quadWeights(q);
                    
                    if EXISTEXACTSOLGRAD
                        
                        % H(div) norm exact
                        norms(f).HDivExact = norms(f).HDivExact + gradUExactVal(comp)^2*detJ*quadWeights(q);
                        
                        % H(div) interpolation norm exact
                        norms(f).HDivExactInt = norms(f).HDivExactInt + gradUExactIntVal(comp)^2*detJ*quadWeights(q);
                        
                        % H(div) error
                        error(f).HDiv = error(f).HDiv + (gradUVal(comp) - gradUExactVal(comp))^2*detJ*quadWeights(q);
                        
                        % H(div) interpolation error
                        error(f).HDivInt = error(f).HDivInt + (gradUVal(comp) - gradUExactIntVal(comp))^2*detJ*quadWeights(q);
                    end
                end
                
                if appCtx.field(f).CALCH1NORM
                    
                    % H1 norm
                    for d=1:dim
                        norms(f).H1 = norms(f).H1 + gradUVal(d)^2*detJ*quadWeights(q);
                    end
                    
                    if EXISTEXACTSOLGRAD
                        
                        % H1 norm exact
                        for d=1:dim
                            norms(f).H1Exact = norms(f).H1Exact + gradUExactVal(d)^2*detJ*quadWeights(q);
                        end
                        
                        % H1 interpolation norm exact
                        for d=1:dim
                            norms(f).H1ExactInt = norms(f).H1ExactInt + gradUExactIntVal(d)^2*detJ*quadWeights(q);
                        end
                        
                        % H1 error
                        for d=1:dim
                            error(f).H1 = error(f).H1 + (gradUVal(d)-gradUExactVal(d))^2*detJ*quadWeights(q);
                        end
                        
                        % H1 interpolation error
                        for d=1:dim
                            error(f).H1Int = error(f).H1Int + (gradUVal(d)-gradUExactIntVal(d))^2*detJ*quadWeights(q);
                        end
                    end
                end
            end
        end
    end
end

% relative errors and sqrt of norms
for f=1:numFields    
    if appCtx.field(f).CALCHDIVNORM
        norms(f).HDiv = sqrt(norms(f).HDiv);
        if EXISTEXACTSOLGRAD
            norms(f).HDivExact = sqrt(norms(f).L2Exact + norms(f).HDivExact);
            norms(f).HDivExactInt = sqrt(norms(f).L2ExactInt + norms(f).HDivExactInt);
            error(f).HDiv = sqrt(error(f).L2 + error(f).HDiv)/norms(f).HDivExact;
            error(f).HDivInt = sqrt(error(f).L2Int + error(f).HDivInt)/norms(f).HDivExactInt;
        end
    end
    if appCtx.field(f).CALCH1NORM
        norms(f).H1 = sqrt(norms(f).H1);
        if EXISTEXACTSOLGRAD
            norms(f).H1Exact = sqrt(norms(f).L2Exact + norms(f).H1Exact);
            norms(f).H1ExactInt = sqrt(norms(f).L2ExactInt + norms(f).H1ExactInt);
            error(f).H1 = sqrt(error(f).L2 + error(f).H1)/norms(f).H1Exact;
            error(f).H1Int = sqrt(error(f).L2Int + error(f).H1Int)/norms(f).H1ExactInt;
        end
    end
    norms(f).L2 = sqrt(norms(f).L2);
    if EXISTEXACTSOL
        norms(f).L2Exact = sqrt(norms(f).L2Exact);
        error(f).L2 = sqrt(error(f).L2/norms(f).L2Exact);
        error(f).L2Int = sqrt(error(f).L2Int/norms(f).L2ExactInt);
    end
end
if EXISTEXACTSOL
    varargout{1} = error;
end