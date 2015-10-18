%% TAYLORHOODVELOCITY
% Taylor-Hood velocity vector basis functions on [-1,1]^2
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
function [numComp, numBasisFuncs, refNodes, basisInfo] = RT0Velocity
numComp = 2;
numBasisFuncs = 4;
refNodes = [ -1 -1 ; 1 -1 ; 1 1 ; -1 1]';

% reference Basis Functions (normalized to [-1,1]^2)
basisInfo(1).fun = @(x,y) [0 ; (y-1)/4];
basisInfo(1).geoType = 'edge';
basisInfo(1).geoNum = 1;
basisInfo(1).DOFType = 'normalFlux';

basisInfo(2).fun = @(x,y) [(x+1)/4 ; 0];
basisInfo(2).geoType = 'edge';
basisInfo(2).geoNum = 2;
basisInfo(2).DOFType = 'normalFlux';

basisInfo(3).fun = @(x,y) [0 ; (y+1)/4];
basisInfo(3).geoType = 'edge';
basisInfo(3).geoNum = 3;
basisInfo(3).DOFType = 'normalFlux';

basisInfo(4).fun = @(x,y) [(x-1)/4 ; 0];
basisInfo(4).geoType = 'edge';
basisInfo(4).geoNum = 4;
basisInfo(4).DOFType = 'normalFlux';