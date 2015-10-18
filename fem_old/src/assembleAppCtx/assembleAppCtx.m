%% ASSEMBLE APPCTX
%
% assemble problem information into appCtx variable
% appCtx =
% field: information specific to a certain variable field
%     name: variable name
%     numComp: number of components of field: 1 for scalars, dim for
%              vectors
%     basisType: type of bais functions used to approximate variable
%     uExact: exact solution, if exists
%     gradUExact: gradient of exact solution, if exists
%     JacForm: Jacobian form
%     bndry: boundary information
%       loc:
%       alpha
%       beta
%       eta
%     nullSpace
%     numBasisFuncs
%     basisInfo
%     numNodeDOF
%     numEdgeDOF
%     numFaceDOF
%     basis
%     basisDer
%     edgeBasis
%     startDOF
%     endDOF
%     cellDOF
%     varForm
%     edgeBC
% dim
% mesh
%     nodes
%     numNodes
%     edges
%     numEdges
%     cells
%     numCells
%     numLocalEdges
%     edgeCells
%     edgeNormals
%     bndryEdges
%     numBndryEdges
% numFields
% maxNumComp
% numCompTotal
% quad
%     prec
%     refNodes
%     refArea
%     v_0Ref
%     numQuadPoints
%     quadPoints
%     quadWeights
%     numEdgeQuadPoints
%     edgeQuadPoints
%     edgeQuadWeights
% numBasisFuncsTotal
% cellGeometry
%     v_0
%     J
%     invJ
%     detJ
%     h
% NewtonRelTol
% NewtonAbsTol
% NewtonMaxIter
% debug
%     displayMesh
%     displayInitialGuess
%     outputResidualError
%     displayResidualError
%     displayCurrentNewtonIterate
%     displayComputedSolution
%     displayExactSolution
%     displayStreamLines
%     displayJacobianSpyPlot
% globalSize

%% VARIABLES---------------------------------------------------------------
if exist('field','var') 
    appCtx.field = field;
    appCtx = assembleVariables(appCtx);
else
    error('myApp:argChk', 'no variables defined')
end

%% DIMENSION---------------------------------------------------------------
appCtx.dim = 2;

%% MESH--------------------------------------------------------------------
if ~exist('x','var')
    if exist('nx','var') && exist('Lx','var')
        x = Lx(1):(Lx(2)-Lx(1))/nx:Lx(2);
    else
        error('myApp:argChk', 'must define either array x or nx and Lx')
    end
end
if ~exist('y','var')
    if exist('ny','var') && exist('Ly','var')
        y = Ly(1):(Ly(2)-Ly(1))/ny:Ly(2);
    else
        error('myApp:argChk', 'must define either array y or ny and Ly')
    end
end
if exist('loc','var')
    appCtx.mesh = assembleMesh(x, y, loc);
else
    appCtx.mesh = assembleMesh(x, y, @(x,y)1); 
end

%% QUADRATURE--------------------------------------------------------------
if exist('prec', 'var')
    appCtx.quad = assembleQuad(prec, appCtx);
else
    error('myApp:argChk', 'quadrature precision prec not defined')
end

%% EXACT SOLUTION----------------------------------------------------------
appCtx = assembleExactSolution(appCtx);

%% FINITE ELEMENT BASIS----------------------------------------------------
appCtx = assembleFE(appCtx);

%% CELL GEOMETRY-----------------------------------------------------------
appCtx.cellGeometry = assembleCellGeometry(appCtx);

%% LOCAL TO GLOBAL MAP-----------------------------------------------------
appCtx = assembleLocalToGlobalMap(appCtx);

%% GLOBAL SIZE-------------------------------------------------------------
appCtx.globalSize = findGlobalSize(appCtx);

%% VARIATIONAL FORM--------------------------------------------------------
if isfield(field, 'varForm')
    appCtx = assembleVariationalForm(appCtx);
else
    error('myApp:argChk', 'variational form not defined')
end

%% JACOBIAN FORM-----------------------------------------------------------
% check whether to use user provided Jacobian form or use a finite
% difference scheme to approximate the Jacobian with parameter h
if ~exist('h','var') || h==0
    if isfield(field,'JacForm')
        appCtx = assembleJacobianForm(appCtx);
        appCtx.h = 0;
	if appCtx.DEBUGLEVEL > 3
	    checkJacobianForm(appCtx);
	end
    else
        error('myApp:argChk', 'Jacobian form must be defined for h ~= 0')
    end
else
    appCtx.h = h;
end

%% BOUNDARY CONDITION------------------------------------------------------
appCtx = assembleBC(appCtx);

%% NEWTON SOLVER-----------------------------------------------------------
if exist('relTol', 'var') && exist('absTol', 'var')
    if exist('maxIter', 'var') 
        appCtx = assembleNewton(relTol,absTol,maxIter,appCtx);
    else
        appCtx = assembleNewton(relTol,absTol,10,appCtx);
        warning('\nno maximum number of Newton iterations provided. maxIter set to 10.\n')
    end
else
    error('myApp:argChk', 'absTol and relTol must be defined for Newton solver')
end

%% NULLSPACE---------------------------------------------------------------
if ~isfield(field, 'nullSpace')
    for f=1:appCtx.numFields
        appCtx.field(f).NULLSPACE = 0;
    end
else 
    appCtx = assembleNullSpace(appCtx);
end

%% STREAM FUNCTION---------------------------------------------------------
if exist('existStreamFun','var') && strcmp(existStreamFun,'true')
    appCtx.EXISTSTREAMFUN = 1;
else
    appCtx.EXISTSTREAMFUN = 0;
end

%% ASSEMBLE FOR C-MEX------------------------------------------------------
% appCtxC = assembleAppCtxC(appCtx);

%% DEGENERATE--------------------------------------------------------------

if exist('degenerate','var') && strcmp(degenerate, 'true')
    appCtx.DEGENERATE = 1;
    if exist('phi','var') && exist('RT0VelocityField','var') && exist('RT0PressureField','var')
        if ~exist('underIntegrate', 'var')
            underIntegrate = 'false';
        end
        if ~exist('Theta','var')
            Theta = 0;
            warning('power Theta unset for degenerate term.  Setting to 0.')
        end
        appCtx = assembleDegenerate(phi, Theta, RT0VelocityField, RT0PressureField ,underIntegrate, appCtx);
    else
        error('myApp:argChk', 'for degenerate problem must define phi, RT0VelocityField, and RT0PressureField')
    end
else
    appCtx.DEGENERATE = 0;
end
appCtx.FORMULATION = FORMULATION;
appCtx.Theta = Theta;
appCtx.phi = phi;
appCtx.dxPhi = dxPhi;

%% CLEAR VARIABLES---------------------------------------------------------
clearvars -except appCtx appCtxC
