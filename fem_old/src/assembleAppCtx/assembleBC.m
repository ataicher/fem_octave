%% ASSEMBLE BC-------------------------------------------------------------
% assemble boundary condition information.
% for each boundary condition on a field f, the user must provide it's
% location on the boundary using a function handle with loc{f}(x,y) = 1
% on the boundary.
% if field f contains a scalar variable, the user provides scalars alpha
% and beta as well as a function handle eta for the robin condition:
%
%   alpha*(natural) + beta*u = eta,
%
% if field f contains a vector variable, the user provides  vectors alpha
% and beta as well as two function handles eta{1} and eta{2} for the robin
% conditions:
%
%   alpha(1)*(natural_normal) + beta(1)*dot(u,n) = eta{1},
%   alpha(2)*(natural_tangent) + beta(2)*dot(u,tau) = eta{2},
%
% where n and tau are the the normal and tangent unit vectors to the
% boundary, respectively. u represent either the scalar or vector
% variables. natural represents the boundary term in the variational form
% multilied by the test function or its normal and tangent parts
%
% INPUT:
%   bndry:  cell array containing the location of each boundary condition.
%           if there are numBC boundary conditions for field f, loc{f}{bc}
%           contains the location of boundary condition bc for field f
%   alpha:  defined above
%   beta:   defined above
%   eta:    defined above
%   appCtx: struct containg problem information
%
% OUTPUT: BC struct containing
%   edgeBCType: numBndryEdges length array containing the boundary
%               condition type for each boundary edge of the mesh
%   bndry:      same as input
%   alpha:      same as input
%   beta:       same as input
%   eta:        same as input
%   numBC:      number of boundary conditions on each field
%
function appCtx = assembleBC(appCtx)

numFields = appCtx.numFields;
numBndryEdges = appCtx.mesh.numBndryEdges;
bndryEdges = appCtx.mesh.bndryEdges;
edges = appCtx.mesh.edges;
nodes = appCtx.mesh.nodes;

if ~isfield(appCtx.field, 'bndry')
    error('myApp:argChk', 'boundary condition not defined')
end

% assign fields that have a boundary condition
for f=1:numFields
    numComp = appCtx.field(f).numComp;
    
    % check which fields have a boundary condition and collect into BC
    % struct
    if ~isempty(appCtx.field(f).bndry)
        appCtx.field(f).BNDRYCONDTION = 1;
        bndry = appCtx.field(f).bndry;
        
        % create edgeBCType array
        numBC = length(bndry);
        edgeBC = constructEdgeBCType(bndry, numBC, numBndryEdges, bndryEdges, nodes, edges);
        for bEPtr=1:numBndryEdges
            if edgeBC(bEPtr) == 0
                plotMesh(appCtx, bndryEdges(bEPtr));
                error('myApp:argChk', 'edge %d does not have a boundary condition', bndryEdges(bEPtr))
            end
            
        end
        
        
        
        for bc = 1:numBC
            
            % check there are enough boundary conditions
            if length(bndry(bc).alpha) ~= numComp || length(bndry(bc).beta) ~= numComp || length(bndry(bc).eta) ~= numComp
                error('myApp:argChk', 'incompatible boundary condition on boundary %d', bc)
            end
            
            % check there are no degenerate boundary conditions
            for comp=1:numComp
                if bndry(bc).alpha(comp) == 0 && bndry(bc).beta(comp) == 0
                    error('myApp:argChk', 'alpha(%d) and beta(%d) cannot both be zero on boundary %d', comp, comp, bc)
                end
            end
            
            % collect into BC struct
            appCtx.field(f).bndry(bc).loc = expandFun(bndry(bc).loc);
            for comp = 1:numComp
                appCtx.field(f).bndry(bc).eta{comp} = expandFun(bndry(bc).eta{comp});
            end
            
        end
        
        % collect into BC struct
        appCtx.field(f).edgeBC = edgeBC;
        appCtx.field(f).numBC = numBC;
        
    else
        appCtx.field(f).BNDRYCONDTION = 0;
    end
    
end
end

%% CONSTRUCTEDGEBCTYPE
% construct edgeBCType array
%
% INPUT:
%   bndry:         cell array containing the location of each boundary
%                  condition. if there are numBC boundary conditions for
%                  field f, loc{f}{bc} contains the location of boundary
%                  condition bc for field f
%   numBC:         number of boundary conditions on for the field
%   numBndryEdges: the number of boundary edges
%   bndyEdges:     numBndryEdges long array. The ith entry is the index of
%                  in edges of the ith boundary edge
%   nodes:         numNodes by 2 array. The ith row contains the (x,y)
%                  coordinates of node_i
%   edges:         numEdges by 2 array. The ith row contains pointers into
%                  nodes, identifying the endpoints of edge_i
%
% OUTPUT:
%   edgeBCType:    numBndryEdges length array containing the boundary
%                  condition type for each boundary edge of the mesh
%
function edgeBC = constructEdgeBCType(bndry, numBC, numBndryEdges, bndryEdges, nodes, edges)

edgeBC = zeros(numBndryEdges,1);
for bEPtr=1:numBndryEdges

    ePtr = bndryEdges(bEPtr);
    n1Ptr = edges(1,ePtr);
    n2Ptr = edges(2,ePtr);
    n1 = nodes(:,n1Ptr);
    n2 = nodes(:,n2Ptr);
    for bc = 1:numBC
        loc = bndry(bc).loc;
        
        if (abs(loc(n1(1),n1(2))) < 1e-10) && (abs(loc(n2(1),n2(2))) < 1e-10)
            edgeBC(bEPtr) = bc;
        end
    end
end
end
