function propTable = expandProp(cfg,dim)
% EXPANDPROP expand the properties table of a config object, so that each
% existing tissue type becomes 'dims' times as many tissue types.
%   Detailed explanation goes here
pre_dims =size(cfg.prop);
assert(pre_dims(1)*dim<3e5,"Requested number of tissue types is larger than 3e5!");
propTable = zeros(pre_dims(1)*dim,pre_dims(2));
for i=1:pre_dims(1)
    propTable((i-1)*dim+1:i*dim,:)=ones(dim,1)*cfg.prop(i,:);
end
end