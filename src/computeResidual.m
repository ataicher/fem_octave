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
function r = computeResidual(u,appCtx)

dim = appCtx.dim;
numCells = appCtx.mesh.numCells;
refArea = appCtx.quad.refArea;
numQuadPoints = appCtx.quad.numQuadPoints;
quadWeights = appCtx.quad.quadWeights';
EXISTEXACTSOL = appCtx.EXISTEXACTSOL;
numFields = appCtx.numFields;
edgeCells = appCtx.mesh.edgeCells;
cells = appCtx.mesh.cells;
edgeNormals = appCtx.mesh.edgeNormals;
numBndryEdges = appCtx.mesh.numBndryEdges;
bndryEdges = appCtx.mesh.bndryEdges;
edges = appCtx.mesh.edges;
nodes = appCtx.mesh.nodes;
numEdgeQuadPoints = appCtx.quad.numEdgeQuadPoints;
globalSize = appCtx.globalSize;
DEGENERATE = appCtx.DEGENERATE;
if DEGENERATE
    phi = appCtx.degenerate.phi;
    RT0VelocityField = appCtx.degenerate.RT0VelocityField;
    RT0PressureField = appCtx.degenerate.RT0PressureField;
    UNDERINTEGRATE = appCtx.degenerate.UNDERINTEGRATE;
    Theta = appCtx.degenerate.Theta;
end

phi = inline('.04 *(x>0) + .01*(x<0)','x','y');
Theta = 0;

% initialize residual vector
r = zeros(globalSize,1);

%% COMPUTE RESIDUAL
for c=1:numCells
    
    % compute transformation from cell c to the reference element
    detJ = appCtx.cellGeometry.detJ{c};
    invJ = appCtx.cellGeometry.invJ{c};
    realArea = detJ*refArea;
    
    % transform quadrature points on the reference element to
    % quadrature points on cell c
    realQuadPoints = projectQuadPoints(c,appCtx);
    
    % push-forward basis functions
    realBasis = projectBasis(c,appCtx);
    realBasisDer = projectBasisDer(c,appCtx);
    %dbstop;
    % get values of u and gradU at the quadrature points of cell c
    [uVal, gradUVal, v_f, gradV_f, p_f, gradP_f, v_m, gradV_m, p_m,gradP_m] = projectDOFLocal(u,c,realBasis,realBasisDer,appCtx);
    %dbstop;
    for f=1:numFields
      numComp = appCtx.field(f).numComp;
      numBasisFuncs = appCtx.field(f).numBasisFuncs;
      rBasis = realBasis.field{f};
      rBasisDer = realBasisDer.field{f};
      cellDOF = appCtx.field(f).cellDOF(:,c);
      varForm = appCtx.field(f).varForm;
      
      % initialize element residual
      elemVec = zeros(numBasisFuncs,1);
      phiValVec = zeros(1,numQuadPoints);
      for q=1:numQuadPoints
	phiValVec(q) = phi(realQuadPoints(1,q),realQuadPoints(2,q));
      end
      	
	if f == 1
	  rBasis_x = reshape(rBasis(:,1,:),numBasisFuncs,numQuadPoints);
	  rBasis_y = reshape(rBasis(:,2,:),numBasisFuncs,numQuadPoints);
	  rBasisDer_xx = reshape(rBasisDer(:,1,1,:),numBasisFuncs,numQuadPoints);
	  rBasisDer_yy = reshape(rBasisDer(:,2,2,:),numBasisFuncs,numQuadPoints);
	  v_f_x = v_f(1,:);
	  v_f_y = v_f(2,:);
	  gradV_f_xx = reshape(gradV_f(1,1,:),1,numQuadPoints);
	  gradV_f_yy = reshape(gradV_f(2,2,:),1,numQuadPoints);
	  
	  elemVec = elemVec + detJ*sum(bsxfun(@times, rBasis_x, quadWeights .* phiValVec.^(-2-2*Theta) .* v_f_x), 2);
	  elemVec = elemVec + detJ*sum(bsxfun(@times, rBasis_y, quadWeights .* phiValVec.^(-2-2*Theta) .* v_f_y), 2);
	  elemVec = elemVec + detJ*sum(bsxfun(@times, rBasisDer_xx, -quadWeights .* p_f), 2);
	  elemVec = elemVec + detJ*sum(bsxfun(@times, rBasisDer_yy, -quadWeights .* p_f), 2);	      

	elseif f == 2
	  rBasis2 = reshape(rBasis(:,1,:), 1, numQuadPoints);
	  
	  elemVec = elemVec + detJ*sum(-quadWeights .* (gradV_f_xx + gradV_f_yy + phiValVec .* (p_f - p_m)));
	  	  
	elseif f == 3
	  rBasis_x = reshape(rBasis(:,1,:), numBasisFuncs, numQuadPoints);
	  rBasis_y = reshape(rBasis(:,2,:), numBasisFuncs, numQuadPoints);
	  rBasisDer_xx = reshape(rBasisDer(:,1,1,:), numBasisFuncs, numQuadPoints);
	  rBasisDer_yy = reshape(rBasisDer(:,2,2,:), numBasisFuncs, numQuadPoints);
	  rBasisDer_xy = reshape(rBasisDer(:,1,2,:), numBasisFuncs, numQuadPoints);
	  rBasisDer_yx = reshape(rBasisDer(:,2,1,:), numBasisFuncs, numQuadPoints);
	  gradV_m_xx = reshape(gradV_m(1,1,:), 1, numQuadPoints);
	  gradV_m_yy = reshape(gradV_m(2,2,:), 1, numQuadPoints);
	  gradV_m_xy = reshape(gradV_m(1,2,:), 1, numQuadPoints);
	  gradV_m_yx = reshape(gradV_m(2,1,:), 1, numQuadPoints);
	
	  elemVec = elemVec + detJ*sum(bsxfun(@times, rBasis_x, -quadWeights .* (1-phiValVec)), 2);
	  elemVec = elemVec + detJ*sum(bsxfun(@times,rBasisDer_xx, quadWeights .* (2*(1-phiValVec) .* gradV_m_xx - (5- 2 * phiValVec)/3 .* (gradV_m_xx + gradV_m_yy) - p_m)),2);		  
	  elemVec = elemVec + detJ*sum(bsxfun(@times,rBasisDer_yy, quadWeights .* (2*(1-phiValVec) .* gradV_m_yy - (5- 2 * phiValVec)/3 .* (gradV_m_xx + gradV_m_yy) - p_m)),2);
	  elemVec = elemVec + detJ * sum(bsxfun(@times,rBasisDer_xy, quadWeights .* ((1-phiValVec) .* (gradV_m_xy + gradV_m_yx))),2);
	  elemVec = elemVec + detJ * sum(bsxfun(@times,rBasisDer_yx, quadWeights .* ((1-phiValVec) .* (gradV_m_xy + gradV_m_yx))),2);
      
	elseif f == 4
	  rBasis4 = reshape(rBasis(:,1,:), 1, numQuadPoints);
	  
	  elemVec = elemVec + detJ*sum(quadWeights.*(-gradV_m_xx - gradV_m_yy + phiValVec .* (p_f - p_m)));
	  
	end   
      
      % insert local residual into global residual vector
      r(cellDOF) = r(cellDOF) + elemVec;
      
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
        
        % <p , phi^(-1/2)*div(phi^(n/2)*v_e)> =
        % (p_c*w_c)*phi_c^(-1/2)*|e|*phi_e^(n/2)*v_e.*n_e
        r(cellDOFV) = r(cellDOFV) - (u(cellDOFP)*rBasisP) * (1/sqrt(phi.c(cDOF))) * eLength .* phi.e(eDOF).^(1+Theta) .* (normalBasisV);
        
        % <phi^(-1/2)*div(phi^(n/2)*u) , w> =
        % w_c*phi_c^(-1/2)* sum_e |e|*phi_e^(n/2)*(u_e*v_e).*n_e)
        r(cellDOFP) = r(cellDOFP) + rBasisP * (1/sqrt(phi.c(cDOF))) * sum( eLength .* phi.e(eDOF).^(1+Theta) .* (u(cellDOFV).*normalBasisV) );
      end
      
      if UNDERINTEGRATE
        %% UNDER-INTEGRATION---------------------------------------------------
        if UNDERINTEGRATE
          cellDOFV = appCtx.field(RT0VelocityField).cellDOF(:,c);
          numCompV = appCtx.field(RT0VelocityField).numComp;
          
          % get value of velocity basis function on each edge
          rBasisV = sign(sum(realBasis.field{RT0VelocityField},3));
          rBasisV = (invJ*(rBasisV'))'/2;
          
          % (u_e,v_e) = (1/2)*|c|*u_e*(v_e.*v_e)
          for comp = 1:numCompV
            r(cellDOFV) = r(cellDOFV) + (1/2) * realArea * (u(cellDOFV).*rBasisV(:,comp)) .* rBasisV(:,comp);
          end
        end
      end
    end
end

%% ELIMINATE CONSTANT NULLSPACE--------------------------------------------
for f=1:numFields
  if appCtx.field(f).NULLSPACE
    nullDOF = appCtx.field(f).NULLSPACE;
    if EXISTEXACTSOL
      uExactSol = projectExactSolution(appCtx);
      r(nullDOF) = u(nullDOF)-uExactSol(nullDOF);
    else
      r(nullDOF) = u(nullDOF);
    end
  end
end

%% BOUNDARY CONDITION
for f=1:numFields
    
  % check if there is a boundary condition on the field
  if appCtx.field(f).BNDRYCONDTION
    numComp = appCtx.field(f).numComp;
    numBasisFuncs = appCtx.field(f).numBasisFuncs;
    edgeBC = appCtx.field(f).edgeBC;
    
    % initialize array of dof eliminated by essential condition
    essentialInd = zeros(globalSize,1); cnt = 1;
    
    for bEPtr=1:numBndryEdges
        
      % determine which boundary condition holds
      bndry = appCtx.field(f).bndry(edgeBC(bEPtr));
      
      % get boundary condition info
      alpha = bndry.alpha;
      beta = bndry.beta;
      eta = bndry.eta;
      
      % determine the cell and local edge
      ePtr = bndryEdges(bEPtr);
      c = edgeCells(1,ePtr);
      lEPtr = find(abs(cells(:,c))==ePtr);
      
      % global dof of cell
      cellDOF = appCtx.field(f).cellDOF(:,c);
      
      % project quadrature and basis functions to real edge
      realEdgeQuadPoints = projectEdgeQuadPoints(lEPtr,c,appCtx);
      realEdgeBasis = projectEdgeBasis(lEPtr,c,appCtx);
      rEdgeBasis = realEdgeBasis.field{f};
      edgeQuadWeights = appCtx.quad.edgeQuadWeights(:,lEPtr);
      
      % determinant for computing edge integrals
      J = appCtx.cellGeometry.J{c};
      e = nodes(:,edges(2,ePtr)) - nodes(:,edges(1,ePtr));
      detJEdge = norm(J*e)/norm(e);
      
      % scalar field
      if numComp == 1
        rEdgeBasis = reshape(rEdgeBasis(:,1,:),numBasisFuncs,numEdgeQuadPoints);
        
        % essential condition
        if alpha == 0
           
          % find nonzero basis function on the boundary
          bInd = sum(abs(rEdgeBasis),2) > 1e-10;
          
          % eliminate basis function from test space
          cntTmp = cnt + nnz(bInd) - 1;
          essentialInd(cnt:cntTmp) = cellDOF(bInd);
          cnt = cntTmp + 1;
          
        % natural condition
        elseif alpha ~= 0 && beta == 0
               
          % evaluate edge integral < eta/alpha , v >
          edgeInt = zeros(numBasisFuncs,1);
          for q=1:numEdgeQuadPoints
            etaVal = eta{1}(realEdgeQuadPoints(1,q),realEdgeQuadPoints(2,q));
            edgeInt = edgeInt + ...
                      detJEdge*edgeQuadWeights(q)*etaVal*rEdgeBasis(:,q)/alpha;
          end
          
          % insert local integrals into global residual
          r(cellDOF) = r(cellDOF) + edgeInt;
          
        % Robin condition
        else
            
          % evaluate u at edge quadrature points
          uLocal = u(cellDOF);
          uVal = uVal + rEdgeBasis'*uLocal;
          
          
          % evaluate edge integral < (eta-beta*u)/alpha , v >
          edgeInt = zeros(numBasisFuncs,1);
          for q=1:numEdgeQuadPoints
            etaVal = eta(realEdgeQuadPoints(1,q),realEdgeQuadPoints(2,q));
            edgeInt = edgeInt + ...
		      detJEdge*edgeQuadWeights(q)*(etaVal - beta*uVal(q))*rEdgeBasis(:,q)/alpha;
          end
          
          % insert local integrals into global residual
          r(cellDOF) = r(cellDOF) + edgeInt;
        end
        
      % vector field
      else
          
        % get normal and tangent direction on edge
        orient = sign(cells(lEPtr,c));
        normal = orient*edgeNormals(:,ePtr);
        tangent = [0 -1 ; 1 0]*normal;
        
        % evalute normal and tangential componenet of basis
        % function at edge quadrature points
        normalBasis = zeros(numBasisFuncs,numEdgeQuadPoints);
        tangentBasis = zeros(numBasisFuncs,numEdgeQuadPoints);
        for q=1:numEdgeQuadPoints
          normalBasis(:,q) = normalBasis(:,q) + rEdgeBasis(:,:,q)*normal;
          tangentBasis(:,q) = tangentBasis(:,q) + rEdgeBasis(:,:,q)*tangent;
        end
        
        % essential condition 1
        if alpha(1) == 0
           
          % find nonzero basis function on the boundary
          bInd = sum(abs(normalBasis),2) > 1e-10;
          
          % eliminate basis function from test space
          cntTmp = cnt + nnz(bInd) - 1;
          essentialInd(cnt:cntTmp) = cellDOF(bInd);
          cnt = cntTmp + 1;
          
        % natural condition 1
        elseif alpha(1) ~= 0 && beta(1) == 0
               
          % evaluate edge integral < eta_1/alpha_1 , v.*n >
          edgeInt = zeros(numBasisFuncs,1);
                    for q=1:numEdgeQuadPoints
                      etaVal = eta{1}(realEdgeQuadPoints(1,q),realEdgeQuadPoints(2,q));
                        edgeInt = edgeInt + ...
                            detJEdge*edgeQuadWeights(q)*etaVal*normalBasis(:,q)/alpha(1);
                    end
                    
                    % insert local integrals into global residual
                    r(cellDOF) = r(cellDOF) + edgeInt;
                    
                % Robin condition 1
                else
                    
                    % evaluate normal componenet of u at edge quadrature points
                    uLocal = u(cellDOF);
                    uNormalVal = normalBasis'*uLocal;
                    
                    % evaluate edge integral < (eta_1-beta_1*u.*n)/alpha_1 , v.*n >
                    edgeInt = zeros(numBasisFuncs,1);
                    for q=1:numEdgeQuadPoints
                        etaVal = eta{1}(realEdgeQuadPoints(1,q),realEdgeQuadPoints(2,q));
                        edgeInt = edgeInt + ...
                            detJEdge*edgeQuadWeights(q)*(etaVal - beta(1)*uNormalVal(q))*normalBasis(:,q)/alpha(1);
                    end
                    
                    % insert local integrals into global residual
                    r(cellDOF) = r(cellDOF) + edgeInt;
                end
                
                % essential condition 2
                if alpha(2) == 0
                        
                    % find nonzero basis function on the boundary
                    bInd = sum(abs(tangentBasis),2) > 1e-10;
                    
                    % eliminate basis function from test space
                    cntTmp = cnt + nnz(bInd) - 1;
                    essentialInd(cnt:cntTmp) = cellDOF(bInd);
                    cnt = cntTmp + 1;
                    
                % natural condition 2
                elseif alpha(2) ~= 0  && beta(2) == 0
                    
                    % evaluate edge integral < eta_2/alpha_2 , v.*tau >
                    edgeInt = zeros(numBasisFuncs,1);
                    for q=1:numEdgeQuadPoints
                        etaVal = eta{2}(realEdgeQuadPoints(1,q),realEdgeQuadPoints(2,q));
                        edgeInt = edgeInt + ...
                            detJEdge*edgeQuadWeights(q)*etaVal*tangentBasis(:,q)/alpha(2);
                    end
                    
                    % insert local integrals into global residual
                    r(cellDOF) = r(cellDOF) + edgeInt;
                    
                % Robin condition 2
                else
                    
                    % evaluate tangential componenet of u at edge quadrature points
                    uLocal = u(cellDOF);
                    uTangentVal = tangentBasis(b,q)'*uLocal;
                    
                    % evaluate edge integral < (eta_2-beta_2*u.*tau)/alpha_2 , v.*tau >
                    edgeInt = zeros(numBasisFuncs,1);
                    for q=1:numEdgeQuadPoints
                        etaVal = eta{2}(realEdgeQuadPoints(1,q),realEdgeQuadPoints(2,q));
                        edgeInt = edgeInt + ...
                            detJEdge*edgeQuadWeights(q)*(etaVal - beta(2)*uTangentVal(q))*tangentBasis(:,q)/alpha(2);
                    end
                    
                    % insert local integrals into global residual
                    r(cellDOF) = r(cellDOF) + edgeInt;
                end
            end
        end
        
        % set residual for essential dof to zero
        r(nonzeros(unique(essentialInd))) = 0;
        
    end
end
