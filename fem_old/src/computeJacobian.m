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
function Jac = computeJacobian(u,appCtx)

dim = appCtx.dim;
numCells = appCtx.mesh.numCells;
refArea = appCtx.quad.refArea;
numQuadPoints = appCtx.quad.numQuadPoints;
quadWeights = appCtx.quad.quadWeights;
h = appCtx.h;
numFields = appCtx.numFields;
edgeCells = appCtx.mesh.edgeCells;
edgeNormals = appCtx.mesh.edgeNormals;
edges = appCtx.mesh.edges;
cells = appCtx.mesh.cells;
numEdgeQuadPoints = appCtx.quad.numEdgeQuadPoints;
numBndryEdges = appCtx.mesh.numBndryEdges;
bndryEdges = appCtx.mesh.bndryEdges;
nodes = appCtx.mesh.nodes;
globalSize = appCtx.globalSize;
DEGENERATE = appCtx.DEGENERATE;
if DEGENERATE
    phi = appCtx.degenerate.phi;
    RT0VelocityField = appCtx.degenerate.RT0VelocityField;
    RT0PressureField = appCtx.degenerate.RT0PressureField;
    UNDERINTEGRATE = appCtx.degenerate.UNDERINTEGRATE;
    Theta = appCtx.degenerate.Theta;
end

% initialize Jacobian matrix without running out of memory
if globalSize^2 < 1e8
    Jac = zeros(globalSize,globalSize);
else
    Jac = sparse(globalSize,globalSize);
end

for c=1:numCells
    
    % get transformation from cell c to the reference element
    detJ = appCtx.cellGeometry.detJ{c};
    invJ = appCtx.cellGeometry.invJ{c};
    realArea = detJ*refArea;
    
    % transform quadrature points on the reference element to
    % quadrature points on cell c
    realQuadPoints = projectQuadPoints(c,appCtx);
    
    % push-forward basis functions
    realBasis = projectBasis(c,appCtx);
    realBasisDer = projectBasisDer(c,appCtx);
    
    % get values of u and gradU at the quadrature points of cell c
    [uVal, gradUVal] = projectDOFLocal(u,c,realBasis,realBasisDer,appCtx);
    
    %% CALCULATE FINITE DIFFERENCE APPROXIMATION TO JACOBIAN
    if h
        for fJ=1:numFields
            numCompJ = appCtx.field(fJ).numComp;
            numBasisFuncsJ = appCtx.field(fJ).numBasisFuncs;
            rBasisJ = realBasis.field{fJ};
            rBasisDerJ = realBasisDer.field{fJ};
            cellDOFJ = appCtx.field(fJ).cellDOF(:,c);
            varForm = appCtx.field(fJ).varForm;
            
            for fI=1:numFields
                numCompI = appCtx.field(fI).numComp;
                numBasisFuncsI = appCtx.field(fI).numBasisFuncs;
                rBasisI = realBasis.field{fI};
                rBasisDerI = realBasisDer.field{fI};
                cellDOFI = appCtx.field(fI).cellDOF(:,c);
                
                elemMat = zeros(numBasisFuncsJ,numBasisFuncsI);
                for q=1:numQuadPoints
                    
                    % apply variational form to u and gradU at each quadrature
                    % point of cell c
                    varFormVVal =  varForm.v(uVal{q},gradUVal{q},realQuadPoints(1,q),realQuadPoints(2,q));
                    varFormGradVVal = varForm.gradV(uVal{q},gradUVal{q},realQuadPoints(1,q),realQuadPoints(2,q));
                    
                    for bI=1:numBasisFuncsI
                        
                        % get values of u and gradU perturbed by h times each basis
                        % function at the quadrature points of cell c
                        uVal_h = uVal{q};
                        gradUVal_h = gradUVal{q};
                        for compI=1:numCompI
                            uVal_h.field{fI}(compI) = uVal_h.field{fI}(compI) + h*rBasisI(bI,compI,q);
                            for d=1:dim
                                gradUVal_h.field{fI}(compI,d) = gradUVal_h.field{fI}(compI,d) + h*rBasisDerI(bI,d,compI,q);
                            end
                        end
                        
                        % apply variational form to perturbed u and gradU at
                        % each quadrature point
                        varFormVVal_h = varForm.v(uVal_h,gradUVal_h,realQuadPoints(1,q),realQuadPoints(2,q));
                        varFormGradVVal_h = varForm.gradV(uVal_h,gradUVal_h,realQuadPoints(1,q),realQuadPoints(2,q));
                        
                        % evaluate integrals ((f0_h-f0)/h , v) + ((f1_h-f1)/h , gradV)
                        for compJ=1:numCompJ
                            elemMat(:,bI) = elemMat(:,bI) + ...
                                detJ*quadWeights(q)*((varFormVVal_h(compJ)-varFormVVal(compJ))/h)*rBasisJ(:,compJ,q);
                            for d=1:dim
                                elemMat(:,bI) = elemMat(:,bI) + ...
                                    detJ*quadWeights(q)*((varFormGradVVal_h(compJ,d)-varFormGradVVal(compJ,d))/h)*rBasisDerJ(:,d,compJ,q);
                            end
                        end
                        
                    end
                end
                
                % insert local Jacobian into global Jacobian matrix
                Jac(cellDOFJ,cellDOFI) = Jac(cellDOFJ,cellDOFI) + elemMat;
            end
        end
        
    %% CALCULATE JACOBIAN FROM JACOBIAN FORM
    else
        for fI=1:numFields
            numCompI = appCtx.field(fI).numComp;
            numBasisFuncsI = appCtx.field(fI).numBasisFuncs;
            rBasisI = realBasis.field{fI};
            rBasisDerI = realBasisDer.field{fI};
            CellDOFI =appCtx.field(fI).cellDOF(:,c);
            JacForm = appCtx.field(fI).JacForm;
            
            for fJ=1:numFields
                numCompJ = appCtx.field(fJ).numComp;
                numBasisFuncsJ = appCtx.field(fJ).numBasisFuncs;
                rBasisJ = realBasis.field{fJ};
                rBasisDerJ = realBasisDer.field{fJ};
                CellDOFJ = appCtx.field(fJ).cellDOF(:,c);
                VW = JacForm.field(fJ).VW;
                VGradW = JacForm.field(fJ).VGradW;
                gradVW = JacForm.field(fJ).gradVW;
                gradVGradW = JacForm.field(fJ).gradVGradW;
                
                elemMat = zeros(numBasisFuncsJ,numBasisFuncsI);
                for q=1:numQuadPoints
                    
                    % compute Jacobian form coeffiecient at u and gradU at each
                    % quadrature point
                    JacFormVW = VW(uVal{q},gradUVal{q},realQuadPoints(1,q),realQuadPoints(2,q));
                    JacFormVGradW = VGradW(uVal{q},gradUVal{q},realQuadPoints(1,q),realQuadPoints(2,q));
                    JacFormGradVW = gradVW(uVal{q},gradUVal{q},realQuadPoints(1,q),realQuadPoints(2,q));
                    JacFormGradVGradW = gradVGradW(uVal{q},gradUVal{q},realQuadPoints(1,q),realQuadPoints(2,q));
                    
                    %evaluate integrals ( df00(u,gradU)*w , v ) +
                    % ( df01(u,gradU)*w , gradV ) + ( df10(u,gradU)*gradW , v )
                    % + ( df11(u,gradU)*gradW , gradV )
                    for compI=1:numCompI
                        for compJ=1:numCompJ
                            elemMat = elemMat + ...
                                detJ*quadWeights(q)*JacFormVW((compI-1)*numCompJ+compJ)*rBasisJ(:,compJ,q)*rBasisI(:,compI,q)';
                            for d=1:dim
                                elemMat = elemMat + ...
                                    detJ*quadWeights(q)*JacFormVGradW((compI-1)*numCompJ*dim+(compJ-1)*dim+d)*rBasisDerJ(:,d,compJ,q)*rBasisI(:,compI,q)';
                                elemMat = elemMat + ...
                                    detJ*quadWeights(q)*JacFormGradVW((compI-1)*dim*numCompJ+(d-1)*numCompJ+compJ)*rBasisJ(:,compJ,q)*rBasisDerI(:,d,compI,q)';
                                for d2=1:dim
                                    elemMat = elemMat + ...
                                        detJ*quadWeights(q)*JacFormGradVGradW((compI-1)*dim*numCompJ*dim+(d-1)*numCompJ*dim+(compJ-1)*dim+d2)*rBasisDerJ(:,d2,compJ,q)*rBasisDerI(:,d,compI,q)';
                                end
                            end
                        end
                    end
                end
                
                % insert local Jacobian into global Jacobian matrix
                Jac(CellDOFJ,CellDOFI) = Jac(CellDOFJ,CellDOFI) + elemMat;
            end
        end
    end
    
    if DEGENERATE
        %% DEGENERATE TERM-----------------------------------------------------
        % get phi values and check that phi.c is nonzero
        cDOF = phi.cellDOF.c(c);
        eDOF = phi.cellDOF.e(:,c);
        if phi.c(cDOF) > 1e-10
            numCompV = appCtx.field(RT0VelocityField).numComp;
            numBasisFuncsV = appCtx.field(RT0VelocityField).numBasisFuncs;
            
            % get value of velocity basis function on each edge
            rBasisV = sign(sum(realBasis.field{RT0VelocityField},3));
            rBasisV = ([invJ(2,2) 0 ; 0 invJ(1,1)]*(rBasisV'))'/2;
            %             rBasisV = (invJ*(rBasisV'))'/2;
            
            % get value of pressure basis function on cell
            rBasisP = realBasis.field{RT0PressureField};
            rBasisP = rBasisP(1,1,1);
            
            % get edge lengths
            ePtr = abs(cells(:,c));
            e = nodes(:,edges(2,ePtr)) - nodes(:,edges(1,ePtr));
            eLength= sqrt(sum(e.^2))';
            
            % get value of v_e.*n_e
            lNormals = [0 -1 ; 1 0 ; 0 1 ; -1 0];
            normalBasisV = zeros(numBasisFuncsV,1);
            for comp = 1:numCompV
                normalBasisV = normalBasisV + rBasisV(:,comp).*lNormals(:,comp);
            end
            
            % get degrees of freedom
            cellDOFV = appCtx.field(RT0VelocityField).cellDOF(:,c);
            cellDOFP = appCtx.field(RT0PressureField).cellDOF(:,c);
            
            % (w_c , phi^(-1/2)*div(phi^(n/2)*v_e) =
            % w_c*phi_c^(-1/2)*|e|*phi_e^(n/2)*v_e.*n_e
            Jac(cellDOFV,cellDOFP) = Jac(cellDOFV,cellDOFP) - rBasisP * (1/sqrt(phi.c(cDOF))) * eLength .* phi.e(eDOF).^(1+Theta) .* normalBasisV;
            
            % (phi^(-1/2)*div(phi^(n/2)*v_e) , w_c) =
            % w_c*phi_c^(-1/2)* |e|*phi_e^(n/2)*v_e.*n_e
            Jac(cellDOFP,cellDOFV) = Jac(cellDOFP,cellDOFV) + ( rBasisP * (1/sqrt(phi.c(cDOF))) * eLength .* phi.e(eDOF).^(1+Theta) .* normalBasisV )';
            
        end
        
        if UNDERINTEGRATE
            %% UNDER INTEGRATION-----------------------------------------------
            cellDOFV = appCtx.field(RT0VelocityField).cellDOF(:,c);
            numCompV = appCtx.field(RT0VelocityField).numComp;
            
            % get value of velocity basis function on each edge
            rBasisV = sign(sum(realBasis.field{RT0VelocityField},3));
            rBasisV = (invJ*(rBasisV'))'/2;
            
            % (v_ei,v_ej) = (1/2)*|c|*(v_ei.*v_je)*delta_ij
            for comp = 1:numCompV
                Jac(cellDOFV,cellDOFV) = Jac(cellDOFV,cellDOFV) + diag( (1/2) * realArea * rBasisV(:,comp) .* rBasisV(:,comp) );
            end
        end
    end
    
end

%% ELIMINATE CONSTANT NULLSPACE--------------------------------------------
for f=1:numFields
    if appCtx.field(f).NULLSPACE
        nullDOF = appCtx.field(f).NULLSPACE;
        nullSpaceVec = zeros(1,globalSize);
        nullSpaceVec(nullDOF) = 1;
        Jac(nullDOF,:) = nullSpaceVec;
        Jac(:,nullDOF) = nullSpaceVec';
    end
end

%% BOUNDARY CONDITIONS-----------------------------------------------------
for fJ=1:numFields
    
    % check if there is a boundary condition on the field
    if appCtx.field(fJ).BNDRYCONDTION
        numCompJ = appCtx.field(fJ).numComp;
        numBasisFuncsJ = appCtx.field(fJ).numBasisFuncs;
        edgeBC = appCtx.field(fJ).edgeBC;
        
        % initialize array of dof eliminated by essential condition
        essentialInd = zeros(globalSize,1); cnt = 1;
        
        for bEPtr=1:numBndryEdges
            
            % determine which boundary condition holds
            bndry = appCtx.field(fJ).bndry(edgeBC(bEPtr));
            
            % get boundary condition info
            alpha = bndry.alpha;
            beta = bndry.beta;
            
            % determine the cell and local edge
            ePtr = bndryEdges(bEPtr);
            c = edgeCells(1,ePtr);
            lEPtr = find(abs(cells(:,c))==ePtr);
            
            % global dof of cell
            cellDOFJ = appCtx.field(fJ).cellDOF(:,c);
            
            % project quadrature and basis functions to real edge
            realEdgeBasis = projectEdgeBasis(lEPtr,c,appCtx);
            rEdgeBasisJ = realEdgeBasis.field{fJ};
            edgeQuadWeights = appCtx.quad.edgeQuadWeights(:,lEPtr);
            
            % determinant for computing edge integrals
            J = appCtx.cellGeometry.J{c};
            e = nodes(:,edges(2,ePtr)) - nodes(:,edges(1,ePtr));
            detJEdge = norm(J*e)/norm(e);
            
            % scalar field
            if numCompJ == 1
                rEdgeBasisJ = reshape(rEdgeBasisJ(:,1,:),numBasisFuncsJ,numEdgeQuadPoints);
                
                % essential condition
                if alpha == 0
                    
                    % find nonzero basis function on the boundary
                    bInd = sum(abs(rEdgeBasisJ),2) > 1e-10;
                    
                    % eliminate basis function from test space
                    cntTmp = cnt + nnz(bInd) - 1;
                    essentialInd(cnt:cntTmp) = cellDOFJ(bInd);
                    cnt = cntTmp + 1;
                    
                    % Robin condition
                elseif beta ~= 0
                    
                    % evaluate edge integral < beta*w/alpha , v >
                    edgeInt = zeros(numBasisFuncsJ,numBasisFuncsJ);
                    for q=1:numEdgeQuadPoints
                        edgeInt = edgeInt + ...
                            detJEdge*edgeQuadWeights(q)*beta*rEdgeBasisJ(:,q)*rEdgeBasisJ(:,q)'/alpha;
                    end
                    
                    % insert local integrals into global Jacobian
                    Jac(cellDOFJ, cellDOFJ) = Jac(cellDOFJ, cellDOFJ) + edgeInt;
                end
                
                % vector field
            else
                
                % get normal and tangent direction on edge
                orient = sign(cells(lEPtr,c));
                normal = orient*edgeNormals(:,ePtr);
                tangent = [0 -1 ; 1 0]*normal;
                
                % evalute normal and tangential componenet of each basis
                % function at edge quadrature points
                normalBasisJ = zeros(numBasisFuncsJ,numEdgeQuadPoints);
                tangentBasisJ = zeros(numBasisFuncsJ,numEdgeQuadPoints);
                for q=1:numEdgeQuadPoints
                    normalBasisJ(:,q) = normalBasisJ(:,q) + rEdgeBasisJ(:,:,q)*normal;
                    tangentBasisJ(:,q) = tangentBasisJ(:,q) + rEdgeBasisJ(:,:,q)*tangent;
                end
                
                % essential condition 1
                if alpha(1) == 0
                    
                    % find nonzero basis function on the boundary
                    bInd = sum(abs(normalBasisJ),2) > 1e-10;
                    
                    % eliminate basis function from test space
                    cntTmp = cnt + nnz(bInd) - 1;
                    essentialInd(cnt:cntTmp) = cellDOFJ(bInd);
                    cnt = cntTmp + 1;
                    
                    % Robin condition 1
                elseif beta(1) ~= 0
                    
                    % evaluate edge integral < beta_1*w.*n/alpha_1 , v.*n >
                    edgeInt = zeros(numBasisFuncsJ,numBasisFuncsJ);
                    for q=1:numEdgeQuadPoints
                        edgeInt = edgeInt + ...
                            detJEdge*edgeQuadWeights(q)*beta(1)*normalBasisJ(:,q)*normalBasisJ(:,q)'/alpha(1);
                    end
                    
                    % insert local integrals into global Jacobian
                    Jac(cellDOFJ, cellDOFJ) = Jac(cellDOFJ, cellDOFJ) + edgeInt;
                end
                
                % essential condition 2
                if alpha(2) == 0
                    
                    % find nonzero basis function on the boundary
                    bInd = sum(abs(tangentBasisJ),2) > 1e-10;
                    
                    % eliminate basis function from test space
                    cntTmp = cnt + nnz(bInd) - 1;
                    essentialInd(cnt:cntTmp) = cellDOFJ(bInd);
                    cnt = cntTmp + 1;
                    
                    % Robin condition 2
                elseif beta(2) ~= 0
                    
                    % evaluate edge integral < beta_2*w.*tau/alpha_2 , v.*tau >
                    edgeInt = zeros(numBasisFuncsJ,numBasisFuncsJ);
                    for q=1:numEdgeQuadPoints
                        edgeInt = edgeInt + ...
                            detJEdge*edgeQuadWeights(q)*beta(2)*tangentBasisJ(:,q)*tangentBasisJ(:,q)'/alpha(2);
                    end
                    
                    % insert local integrals into global Jacobian
                    Jac(cellDOFJ, cellDOFJ) = Jac(cellDOFJ, cellDOFJ) + edgeInt;
                end
            end
        end
        
        % set Jacobian rows and columns corresponding to essential
        % condition to the identity
        essentialInd = nonzeros(unique(essentialInd))';
        Jac(essentialInd,:) = 0;
        Jac(:,essentialInd) = 0;
        for ind=essentialInd
            Jac(ind,ind) = 1;
        end
        
    end
end