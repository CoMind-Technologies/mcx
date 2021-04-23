function [cfg1,outputArg2] = updateProperties(cfg0,propVolMatrix)
%UPDATEPROPERTIES Update the properties of a config file retrieved from
% 'propVolMartix'. Various things have to match-up:
%       'propVolMatrix' must have the same # columns as cfg0.prop (if it
%       exists!)
%       NB: Currently, this file can ONLY deal with changes in 
if isfield(cfg0,'prop')
    assert(size(cfg0.prop,2)==size(propVolMatrix,2),"Mismatch in column number between cfg0.prop and propVolMatrix.")
end
if isfield(cfg0,'propVol')
   % config has a propVol matrix. Update it. 
end
cfg1 = cfg0;
outputArg2 = inputArg2;
end

