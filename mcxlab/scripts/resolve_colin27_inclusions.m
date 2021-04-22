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
cfg0.nphoton=5e5; % number of photons
cfg0.tstart=0;
cfg0.tend=5e-9;
cfg0.tstep=2e-10;
cfg0.isreflect=1; % enable reflection at exterior boundary. Kinda cheating
cfg0.isrefint=1;  % enable reflection at interior boundary too
cfg0.issavedet=0; % Record partial pathlength of detected photons
cfg0.ismomentum=1; % Record momenta of detected photons
cfg0.vol=colin27; % cfg0 is 'reality': incl. inclusions:
cfg0.srcpos=[149 75 149];
cfg0.srcdir=normalize([-1 0 -1],'norm');
cfg0.detpos=[137 75 160 4
             149 65 145 4];
cfg0.prop=[  0         0         1.0000    1.0000 % (0) background/air
            0.0190    7.8182    0.8900    1.3700 % (1) scalp
            0.0190    7.8182    0.8900    1.3700 % (2) skull
            0.0040    0.0090    0.8900    1.3700 % (3) csf
            0.0200    9.0000    0.8900    1.3700 % (4) gray matters
            0.0800   40.9000    0.8400    1.3700 % (5) white matters
                 0         0    1.0000    1.0000]; % (6) air pockets
% GPU Config
cfg0.autopilot=1;
cfg0.gpuid=1;
% Now create copies; for nwo just 1
cfg_incl=cfg0;
% Add the inclusion and give it a type number
cfg_incl.vol(csf_inclusions(set_indices(1)).indices)=max(max(max(colin27)))+1;
% Ammend the properties spec:
cfg_incl.prop=[cfg0.prop; .01  0.0090    0.8900    1.3700]; % Only pert absorbtion
% mcxpreview(cfg1)
%% Run baseline; i.e. with inclusion
fprintf('Running baseline simulation ...\n');
[fluence,detphotons,~,seeds]=mcxlab(cfg_incl);
%% Now get ABSORPTION jacobians
cfg0.outputtype='jacobian';
cfg0.detphotons=detphotons.data;
cfg0.seed=seeds.data;
fprintf('Computing Jacobian...\n');
[flux2, detp2, vol2, seeds2]=mcxlab(cfg0);
jacobian = sum(flux2,4); % Summing over time gates
%% Now update the baseline assumption...
% Find all non-zero vozels of the jacobian:
effected_tissues=cfg0.vol(find(jacobian));