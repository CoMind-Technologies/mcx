clear all;
clear cfg0;
try
    gpuinfo=mcxlab('gpuinfo');
catch
    USE_MCXCL=1;
end
%% For reproducibility
rng('default')
%% load volume and inclusions
load colin27_v3.mat
load colin27_csf_500vol.mat
%% Cherry pick the ones closest to the crown
maxes= [0 0 0 0];
set_indices = [0 0 0 0];
for i=1:length(csf_inclusions)
    max_z = max(csf_inclusions(i).voxels(:,3));
    if maxes(end)<max_z
        maxes(1:3)=maxes(2:4);
        set_indices(1:3)=set_indices(2:4);
        maxes(end) = max_z;
        set_indices(end)=i;
    end
end
%% Create a basline simulation config file and a set for each inclusion
cfg0.vol=colin27; % cfg0 is 'reality': incl. inclusions:
cfg0.prop=[  0         0         1.0000    1.0000 % (0) background/air
            0.0190    7.8182    0.8900    1.3700 % (1) scalp
            0.0190    7.8182    0.8900    1.3700 % (2) skull
            0.0040    0.0090    0.8900    1.3700 % (3) csf
            0.0200    9.0000    0.8900    1.3700 % (4) gray matters
            0.0800   40.9000    0.8400    1.3700 % (5) white matters
                 0         0    1.0000    1.0000]; % (6) air pockets

% Now create copies for inclusions
cfg_incl=cfg0;
cfg_incl.nphoton=5e5; % number of photons
cfg_incl.srcpos=[149 75 149];
cfg_incl.srcdir=normalize([-1 0 -1],'norm');
cfg_incl.detpos=[137 75 160 4
             149 65 145 4];
cfg_incl.tstart=0;
cfg_incl.tend=5e-9;
cfg_incl.tstep=2e-10;
cfg_incl.isreflect=1; % enable reflection at exterior boundary. Kinda cheating
cfg_incl.isrefint=1;  % enable reflection at interior boundary too
cfg_incl.issavedet=0; % Record partial pathlength of detected photons
cfg_incl.ismomentum=1; % Record momenta of detected photons
% Add the inclusion and give it a type number
cfg_incl.vol(csf_inclusions(set_indices(1)).indices)=max(max(max(colin27)))+1;
% Ammend the properties spec:
cfg_incl.prop=[cfg0.prop; .01  0.0090    0.8900    1.3700]; % Only pert absorbtion
% GPU Config
cfg_incl.autopilot=1;
cfg_incl.gpuid=1;
mcxpreview(cfg_incl)
%% Run baseline; i.e. with inclusion
fprintf('Running baseline simulation ...\n');
[fluence,detphotons,~,seeds]=mcxlab(cfg_incl);
cfg_incl.dephotons=detphotons; % Add detected photons to config file
%% Save both
save('./configs/colin27_cfg.mat','cfg0');
save('./configs/colin27_cfs_inclusion_cfg.mat','cfg_incl');