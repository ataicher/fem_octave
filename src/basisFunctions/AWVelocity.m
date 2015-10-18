%% AWVELOCITY
% The AW_0 basis functions (on an element) are of the following types
%   vxcSS = vel x-dir on corner SS
%   vxbS  = vel x-dir on edge S
% 
%   vycSS = vel y-dir on corner SS
%   vybS  = vel y-dir on edge S
% 
%   where S is one of
%     M = minus (to the left or down)
%     P = plus (to the right or up)
% 
% indices proceed counter-clockwise over nodes beginning in the 
% bottom-left, then counter-clockwise over edges beginning with bottom.
%
%   vxcMM, vycMM, vxcPM, VycPM, vxcMP, vycMP, vxcPP, vycPP,
%   vybM, vxbP, vybP, vxbM
% 
%  vxcMP,vycMP-----vybP-----vxcPP,vycPP
%   |                             |
%   |                             |
%  vxbM                          vxbP  
%   |                             |
%   |                             |
%  vxcMM,vycMM-----vybM-----vxcPM,vycPM
% 
% vx and vy basis functions are vector valued.  The naming for the x-basis
% function and its x-component is overloaded. Similarly for the y-basis 
% functions. i.e.,
% 
%          vxcSS = / vxcSS \   vycSS = /   0   \
%                  \   0   / ,         \ vycSS /
% 
% the bottom left corner and left bubble nodal x-basis functions on [-1,1]^2 are defined by the relations
% 
%   vxcMM(-1,-1) = 1
%   vxcMM(-1,1)  = 0
%   vxcMM(1,-1)  = 0
%   vxcMM(1,1)   = 0
%   average normal component on each edge is 0
% 
%   vxbM(-1,-1) = 0
%   vxbM(-1,1)  = 0
%   vxbM(1,-1)  = 0
%   vxbM(1,1)   = 0
%   average normal component on left edge is 1
%   average normal component on remaining edges is 0
% 
% basis functions are given by
%   vxcMM(x,y) = (1/8)*(1-x)*(y-1)*(1+3y)
%   vxbM(x,y)  = (3/4)*(1-x)*(1-y^2)
%
% the remainder of the basis functions are defined using symmetries.
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
function [numComp, numBasisFuncs, refNodes, basisInfo] = AWVelocity
numComp = 2;
numBasisFuncs = 12;
refNodes = [ -1 -1 ; 1 -1 ; 1 1 ; -1 1]';

% reference Basis Functions (normalized to [-1,1]^2)
basisInfo(1).fun = @(x,y) [(1-x)*(y-1)*(1+3*y)/8 ; 0];
basisInfo(1).gradFun = @(x,y) [(1-y)*(1+3*y)/8 -(1-x)*(1-3*y)/4 ; 0 0]; 
basisInfo(1).geoType = 'node';
basisInfo(1).geoNum = 1;
basisInfo(1).DOFType = 'nodal1';

basisInfo(2).fun = @(x,y) [0 ; (1-y)*(x-1)*(1+3*x)/8];
basisInfo(2).gradFun = @(x,y) [0 0 ; -(1-y)*(1-3*x)/4 -(1-x)*(1-3*x)/8]; 
basisInfo(2).geoType = 'node';
basisInfo(2).geoNum = 1;
basisInfo(2).DOFType = 'nodal2';

basisInfo(3).fun = @(x,y) [(1+x)*(y-1)*(1+3*y)/8 ; 0];
basisInfo(3).gradFun = @(x,y) [(y-1)*(1+3*y)/8 -(1+x)*(1-3*y)/4 ; 0 0]; 
basisInfo(3).geoType = 'node';
basisInfo(3).geoNum = 2;
basisInfo(3).DOFType = 'nodal1';

basisInfo(4).fun = @(x,y) [0 ; (1-y)*(-x-1)*(1-3*x)/8];
basisInfo(4).gradFun = @(x,y) [0 0 ; (1-y)*(1+3*x)/4 (1+x)*(1-3*x)/8]; 
basisInfo(4).geoType = 'node';
basisInfo(4).geoNum = 2;
basisInfo(4).DOFType = 'nodal2';

basisInfo(5).fun = @(x,y) [(1+x)*(-y-1)*(1-3*y)/8; 0];
basisInfo(5).gradFun = @(x,y) [ -(1+y)*(1-3*y)/8 (1+x)*(1+3*y)/4 ; 0 0]; 
basisInfo(5).geoType = 'node';
basisInfo(5).geoNum = 3;
basisInfo(5).DOFType = 'nodal1';

basisInfo(6).fun = @(x,y) [0 ; (1+y)*(-x-1)*(1-3*x)/8];
basisInfo(6).gradFun = @(x,y) [0 0 ; (1+y)*(1+3*x)/4 -(1+x)*(1-3*x)/8]; 
basisInfo(6).geoType = 'node';
basisInfo(6).geoNum = 3;
basisInfo(6).DOFType = 'nodal2';

basisInfo(7).fun = @(x,y) [(1-x)*(-y-1)*(1-3*y)/8 ; 0];
basisInfo(7).gradFun = @(x,y) [(1+y)*(1-3*y)/8 (1-x)*(3*y+1)/4 ; 0 0];
basisInfo(7).geoType = 'node';
basisInfo(7).geoNum = 4;
basisInfo(7).DOFType = 'nodal1';

basisInfo(8).fun = @(x,y) [0; (1+y)*(x-1)*(1+3*x)/8];
basisInfo(8).gradFun = @(x,y) [0 0 ; -(1+y)*(1-3*x)/4 (x-1)*(1+3*x)/4];
basisInfo(8).geoType = 'node';
basisInfo(8).geoNum = 4;
basisInfo(8).DOFType = 'nodal2';

basisInfo(9).fun = @(x,y) [0 ; -3*(1-y)*(1-x^2)/8];
basisInfo(9).gradFun = @(x,y) [0 0 ; 3*(1-y)*x/4 3*(1-x^2)/8];
basisInfo(9).geoType = 'edge';
basisInfo(9).geoNum = 1;
basisInfo(9).DOFType = 'normalFlux';

basisInfo(10).fun = @(x,y) [3*(1+x)*(1-y^2)/8 ; 0];
basisInfo(10).gradFun = @(x,y) [3*(1-y^2)/8 -3*(1+x)*y/4 ; 0 0];
basisInfo(10).geoType = 'edge';
basisInfo(10).geoNum = 2;
basisInfo(10).DOFType = 'normalFlux';

basisInfo(11).fun = @(x,y) [0 ; 3*(1+y)*(1-x^2)/8];
basisInfo(11).gradFun = @(x,y) [0 0 ; -3*(1+y)*x/4 3*(1-x^2)/8];
basisInfo(11).geoType = 'edge';
basisInfo(11).geoNum = 3;
basisInfo(11).DOFType = 'normalFlux';

basisInfo(12).fun = @(x,y) [-3*(1-x)*(1-y^2)/8 ; 0];
basisInfo(12).gradFun = @(x,y) [3*(1-y^2)/8 3*(1-x)*y/4 ; 0 0];
basisInfo(12).geoType = 'edge';
basisInfo(12).geoNum = 4;
basisInfo(12).DOFType = 'normalFlux';
