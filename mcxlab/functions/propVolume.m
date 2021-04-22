function [outputArg1,outputArg2] = propVolume(cfg,varargin)
%PROPVOLUME: return a volume_id x properties matrix
%   The given config file (cfg) must have a .prop matrix and a .vol volume file,
%   as per the 'usual' mcx specifications.
%   Optional input is which property colums to return. Must be saame/less than #
%   of properties (<=4).
% RETURNS
X=product(size(config.vol));
if nargin>0
    outputArg1 = zeros(X,nargin);
else
    outputArg1 = zeros(X,4);
end
for i=1:X
    outputArg1(i,:)=cfg.prop(config.vol(i),:);
end

