%% PIECEWISECONST
% basis function with constant value of 1 over each element
%
% INPUT:
%
% OUTPUT:
%   numComp:       number of components of the basis functions (vector or
%                  scalar)
%   numBasisFuncs: number of basis functions
%   refNodes:      nodes of reference element
%   nodeDOFType:   degree of freedom type on each node (same for each node)
%                  accepts degrees of freedom of the following types: 
%                  nodal: nodal value for scalar-valued basis
%                  nodal1: x-component value for vector-valued basis
%                  nodal2: y-component value for vector-valued basis
%   edgeDOFType:   degree of freedom type on each edge (same for each edge)
%                  accepts degrees of freedom of the following types:
%                  normalFlux: normal flux across edge
%                  tangentialFlux: tangential flux across edge
%   faceDOFType:   degree of freedom type on each face (same for each face)
%                  accepts degrees of freedom of the following types:
%                  averageVal: average value on the face
%   basisFuncs:    numComp by numBasisFuncs cell array of strings with 
%                  equations for each component of each basis function. 
%                  basisFuncs{i,j} contains a string representing the
%                  jth component of the ith basis function. basis funcions 
%                  are ordered by a local numbering convention:
%                  1. nodes are numbered first, then edges, then faces
%                  2. numbering begins with the bottom left dof and 
%                  continues in a counter-clockwise direction
%
function [dim, numComp, numBasisFuncs, refNodes, basisInfo] = piecewiseConst
% 
dim = 1;
numComp = 1;
numBasisFuncs = 2;
refNodes = [ -1 -1]';

% reference Basis Functions (normalized to [-1,1]^2)
basisInfo(1).fun = @(x,y) (1-x)/2;
basisInfo(1).geoType = 'edge';
basisInfo(1).geoNum = 1;
basisInfo(1).DOFType = 'normal';

basisInfo(1).fun = @(x,y) (x-1)/2;
basisInfo(1).geoType = 'edge';
basisInfo(1).geoNum = 2;
basisInfo(1).DOFType = 'normal';
