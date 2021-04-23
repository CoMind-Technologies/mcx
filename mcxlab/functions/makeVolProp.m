function volPropMatrix = makeVolProp(config,varargin)
%PROPVOLUME: return a volume_id x properties matrix
%   The given config file (cfg) must have a .prop matrix and a .vol volume file,
%   as per the 'usual' mcx specifications.
%   Optional input is which property colums to return. Must be same/less than #
%   of properties (<=4).
% RETURNS
X=prod(size(config.vol));
if length(varargin)>0
    volPropMatrix = zeros(X,length(varargin));
    for i=1:X
        volPropMatrix(i,cell2mat(varargin))=config.prop(config.vol(i)+1,cell2mat(varargin));
    end
else 
    volPropMatrix = zeros(X,4); % Assuming the default 4 properties.
    for i=1:X
        volPropMatrix(i,:)=config.prop(config.vol(i)+1,:);
    end
end