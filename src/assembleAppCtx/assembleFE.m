%% ASSEMBLEFE
% assemble information about the finite element basis
% local numbering convention:
%   1. nodes are numbered first, then edges, then faces
%   2. numbering begins with the bottom left dof and continues in a
%      counter-clockwise direction
%
% INPUT:
%  basisType: function handle to Matlab function that outputs in order:
%             numComp:       number of components of basis functions
%             numBasisFuncs: number of basis functions
%             refNodes:      reference nodes of reference element
%             nodeDOFType:   degree of freedom type on each node
%             edgeDOFType:   degree of freedom type on each edge
%             faceDOFType:   degree of freedom type on each face
%             basisFuncs:    numComp by numBasisFuncs cell array of string
%                            with equations for each component of each
%                            basis function. basisFuncs{i,j} contains
%                            a string representing the ith component of the
%                            jth basis function
%  appCtx:    struct containing problem information
%
% OUTPUT: FE struct containing
%   basisType:     same as input
%   numBasisFuncs: number of basis functions for each field
%   nodeDOFType:   number of degrees of freedom on each node for each field
%   edgeDOFType:   number of degrees of freedom on each edge for each field
%   faceDOFType:   number of degrees of freedom on each face for each field
%   nodeDOFType:   degree of freedom type on each node
%   edgeDOFType:   degree of freedom type on each edge
%   faceDOFType:   degree of freedom type on each face
%   basisFuncs:    function handles containing the basis functions.
%                  basisfuncs{comp,b} contains component comp of the basis
%                  function b
%   basis:         array containing basis function values at quadrature
%                  points on the reference element for each field.
%                  basis{f}(b,comp,q) contains the value of basis function
%                  b for component comp at quadrature points q of field f
%   basisDer:      array containing basis function derivative values at
%                  quadrature points on the reference element for each
%                  field. basisDer(b,d,comp,q) contains the value of basis
%                  function b derivative in the direction d for component
%                  comp at quadrature points q of field f
%   edgeBasis      array containing basis function values at quadrature
%                  points on the reference element edges for each field.
%                  edgeBasis{f}(b,comp,q,lE) contains the value of basis
%                  function b for component comp at quadrature points q of
%                  local edge lE of field f
%   edgeBasisDer:  array containing basis function derivative values at
%                  quadrature points on the reference element edges for
%                  each field. edgeBasisDer(b,d,comp,q,lE) contains the
%                  value of basis function b derivative in the direction d
%                  for component comp at quadrature points q of local edge
%                  lE of field f
%   v_0Ref:        bottom left vertex coordinates of reference element
%   refArea:       area of reference element
%
function appCtx = assembleFE(appCtx)

numFields = appCtx.numFields;
numLocalEdges = appCtx.mesh.numLocalEdges;
prec = appCtx.quad.prec;

% check basisType is defined
if ~isfield(appCtx.field, 'basisType') 
    error('myApp:argChk', 'finite element basis not defined')
end

L2BasisNames = { 'piecewiseConst'};
HDivBasisNames = {'RT0Velocity'};
H1BasisNames = {'AWVelocity', 'P1Quad', 'P1Quad2','TaylorHoodVelocity'};

numBasisFuncsTotal = 0;
for f=1:numFields
    numComp = appCtx.field(f).numComp;
    
    % check basisType is defined
    if isa(appCtx.field(f).basisType,'function_handle')
        CALCHDIVNORM = 0;
        CALCH1NORM = 0;
        if max(strcmp(func2str(appCtx.field(f).basisType), L2BasisNames))
        elseif max(strcmp(func2str(appCtx.field(f).basisType), HDivBasisNames))
            CALCHDIVNORM = 1;
            
        elseif max(strcmp(func2str(appCtx.field(f).basisType), H1BasisNames))
            CALCH1NORM = 1;
        else
            error('myApp:argChk', 'finite element basis for field %d not defined', f)
        end
    else
        error('myApp:argChk', 'finite element basis for field %d must be a function handle', f)
    end
        
        
    
    % obtain basis function information
    [basisNumComp, numBasisFuncs, ~, basisInfo] = appCtx.field(f).basisType();

    % check consitency between number of components of variables and basis
    % function
    if numComp ~= basisNumComp
        error('myApp:argChk', 'Incompatible number of components for finite element basis of field %d', f)
    end
    
    % number of degrees of freedom for each node, edge, and face
    [numNodeDOF, numEdgeDOF, numFaceDOF] = getDOFNumbers(numBasisFuncs, numLocalEdges, basisInfo);
    
    % load basis, basisDer, and edgeBasis arrays
    load(['./src//basisFunctions/assembledBasis/', func2str(appCtx.field(f).basisType), num2str(prec)]);
    
    % total number of basis functions
    numBasisFuncsTotal = numBasisFuncsTotal + numBasisFuncs;
    
    % collect into FE struct
    appCtx.field(f).basisType = func2str(appCtx.field(f).basisType);
    appCtx.field(f).numBasisFuncs = numBasisFuncs;
    appCtx.field(f).basisInfo = basisInfo;
    appCtx.field(f).numNodeDOF = numNodeDOF;
    appCtx.field(f).numEdgeDOF = numEdgeDOF;
    appCtx.field(f).numFaceDOF = numFaceDOF;
    appCtx.field(f).basis = basis;
    appCtx.field(f).basisDer = basisDer;
    appCtx.field(f).edgeBasis = edgeBasis;
    appCtx.field(f).CALCHDIVNORM = CALCHDIVNORM;
    appCtx.field(f).CALCH1NORM = CALCH1NORM;
end

appCtx.numBasisFuncsTotal = numBasisFuncsTotal;
end

%% GETDOFNUMBERS
function [numNodeDOF, numEdgeDOF, numFaceDOF] = getDOFNumbers(numBasisFuncs, numLocalEdges, basisInfo)

    numNodeDOF = zeros(numLocalEdges,1);
    numEdgeDOF = zeros(numLocalEdges,1);
    numFaceDOF = 0;
    for b=1:numBasisFuncs
        
        % nodes
        if strcmp(basisInfo(b).geoType, 'node')
            lNPtr = basisInfo(b).geoNum;
            numNodeDOF(lNPtr) = numNodeDOF(lNPtr) + 1;
            
        % edges
        elseif strcmp(basisInfo(b).geoType, 'edge')
            lEPtr = basisInfo(b).geoNum;
            numEdgeDOF(lEPtr) = numEdgeDOF(lEPtr) + 1;
            
        %faces
        elseif strcmp(basisInfo(b).geoType, 'face')
            numFaceDOF = numFaceDOF + 1;
        end
    end
    
    % check there are the same number of degrees of freedom on each node
    if length(unique(numNodeDOF)) == 1
        numNodeDOF = numNodeDOF(1);
    else
        error('myApp:argChk', 'number of dof of field %d must be equal for each node', f)
    end
    
    % check there are the same number of degrees of freedom on each edge
    if length(unique(numEdgeDOF)) == 1
        numEdgeDOF = numEdgeDOF(1);
    else
        error('myApp:argChk', 'number of dof of field %d must be equal for each edge', f)
    end
end