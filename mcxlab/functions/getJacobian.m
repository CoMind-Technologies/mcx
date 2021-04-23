function [jacobian,mu_a_vector] = getJacobian(cfgs)
% GETJACOBIANS: Get the Jacobians of these configuration files, 
% base config.
%   'cfgs' is a config file (or list of config files) each with all run
%   properties and with detected photon objects AND seeds!
%
try
    gpuinfo=mcxlab('gpuinfo');
catch
    USE_MCXCL=1;
end
cfgs.outputtype='jacobian';
flux=mcxlabcl(cfgs)
jacobian=sum(flux.data,4);
try
    mu_a_vector=cfgs.volprop(:,1);
catch
    vp = volProp(cfg);
    mu_a_vector=vp(:,1);
end