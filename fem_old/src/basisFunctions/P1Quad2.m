%% P1QUAD2
% P1 nodal vector basis functions on [-1,1]^2
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
function [numComp, numBasisFuncs, refNodes, basisInfo] = P1Quad2

numComp = 2;
numBasisFuncs = 8;
refNodes = [ -1 -1 ; 1 -1 ; 1 1 ; -1 1]';

basisInfo(1).fun = @(x,y) [1/4*(x - 1)*(y - 1) ; 0];
basisInfo(1).geoType = 'node';
basisInfo(1).geoNum = 1;
basisInfo(1).DOFType = 'nodal1';

basisInfo(2).fun = @(x,y) [0 ; 1/4*(x - 1)*(y - 1)];
basisInfo(2).geoType = 'node';
basisInfo(2).geoNum = 1;
basisInfo(2).DOFType = 'nodal2';

basisInfo(3).fun = @(x,y) [1/4*(-x - 1)*(y - 1) ; 0];
basisInfo(3).geoType = 'node';
basisInfo(3).geoNum = 2;
basisInfo(3).DOFType = 'nodal1';

basisInfo(4).fun = @(x,y) [0 ; 1/4*(-x - 1)*(y - 1)];
basisInfo(4).geoType = 'node';
basisInfo(4).geoNum = 2;
basisInfo(4).DOFType = 'nodal2';

basisInfo(5).fun = @(x,y) [1/4*(-x - 1)*(-y - 1) ; 0];
basisInfo(5).geoType = 'node';
basisInfo(5).geoNum = 3;
basisInfo(5).DOFType = 'nodal1';

basisInfo(6).fun = @(x,y) [0 ; 1/4*(-x - 1)*(-y - 1)];
basisInfo(6).geoType = 'node';
basisInfo(6).geoNum = 3;
basisInfo(6).DOFType = 'nodal2';

basisInfo(7).fun = @(x,y) [1/4*(x - 1)*(-y - 1) ; 0];
basisInfo(7).geoType = 'node';
basisInfo(7).geoNum = 4;
basisInfo(7).DOFType = 'nodal1';

basisInfo(8).fun = @(x,y) [0 ; 1/4*(x - 1)*(-y - 1)];
basisInfo(8).geoType = 'node';
basisInfo(8).geoNum = 4;
basisInfo(8).DOFType = 'nodal2';