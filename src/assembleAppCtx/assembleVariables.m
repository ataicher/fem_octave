%% ASSEMBLE VARIABLES------------------------------------------------------
% Assemble variables information
% INPUT:
%   numFields: total number of fields. e.g. velocity, pressure, etc.
%   numComp:   number of components in the field. i.e. numComp{i} is the
%              number of components of field_i
%
% OUTPUT: variables struct containing
%   numFields:    same as input
%   numComp:      same as input
%   naxNumComp:   max number of component for any field
%   numCompTotal: sum of the number of components for all fields
%
function appCtx = assembleVariables(appCtx)

% obtain number of fields
numFields = length(appCtx.field);

% create and resize names cell array
if isfield(appCtx.field, 'name')
    for f =1:numFields
        if isempty(appCtx.field(f).name)
            appCtx.field(f).name = ['field ', num2str(f)];
        end
    end
    
else
    for f=1:numFields
        appCtx.field(f).name = ['field ', num2str(f)];
    end
end

maxNumComp = 0;
numCompTotal = 0;
for f=1:numFields;
    
    % check numComp is defined
    if ~isfield(appCtx.field(f), 'numComp')
        error('myApp:argChk', 'numComp undefined for field %d', f);
    end
    
    numComp = appCtx.field(f).numComp;
    
    % check numComp is 1 or 2
    if ~(numComp == 1 || numComp == 2)
        error('myApp:argChk', 'numComp must be 1 for scalar fields or 2 for vector fields');
    end
    
    % find maxNumComp and numCompTotal
    numCompTotal = numCompTotal + numComp;
    maxNumComp = max(maxNumComp,numComp);
end

%% COLLECT INTO VARIABLES STRUCT-------------------------------------------
appCtx.numFields = numFields;
appCtx.maxNumComp = maxNumComp;
appCtx.numCompTotal = numCompTotal;
end