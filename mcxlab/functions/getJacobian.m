function [Jacobians] = insert_inclusions(cfgs,base_cfg)
% GETJACOBIANS: Get the jacobians between the pre-run config files and the
% base config.
%   'cfgs' is a config file (or list of config files) each with all run
%   properties and with detected photon objects AND seeds!
%   
temp_cfgs = cfgs;
for cfg_=temp_cfgs
    cfg_.vol=base_cfg.vol;
    cfg_.prop=base
end

