%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tdNIRS simulation in colin27 voxel tissue model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clear cfg;
try
    gpuinfo=mcxlab('gpuinfo');
catch
    USE_MCXCL=1;
end

%% Define simulation parameters
    % making a minimal set of simulation parameters
% set seed to make the simulation repeatible
cfg.seed=hex2dec('623F9A9E'); 
% Volume model being used.
load colin27_v3.mat
cfg.vol=colin27;
cfg.nphoton=1e6;
% Each line below defines the optical parameters for that tissue type:
% mua (absorption) mus (scattering) g (anisotropy) and 'n'
cfg.prop=[  0         0         1.0000    1.0000 % (0) background/air
        0.0190    7.8182    0.8900    1.3700 %(1) scalp
        0.0190    7.8182    0.8900    1.3700 %(2) skull
        0.0040    0.0090    0.8900    1.3700 %(3) csf
        0.0200    9.0000    0.8900    1.3700 %(4) gray matters
        0.0800   40.9000    0.8400    1.3700 %(5) white matters
             0         0    1.0000    1.0000]; %(6) air pockets
%
cfg.srcpos=[75 14 73]; %Source position
cfg.detpos=[75 25. 92 10]; %Detector Position    
cfg.srcdir=[0.0 1.0 0.0, 5.0];% define the source DIRECTION
% time-domain simulation parameters
cfg.tstart=0;% Start
cfg.tend=5e-9;% Stop
cfg.tstep=2e-10;% timestep  
%% GPU thread configuration
cfg.autopilot=1;% Let MCX decide how many threads are important
cfg.gpuid=1;% 

%% Run Simulation
fprintf('running simulation ...\n');
[out.fluence,out.detphoton,out.vol,out.seed,out.trajectory] = mcxlab(cfg);
%% Is the mismatch?!
    % The volume structures are 3-dimensional, so let's see if ther are ANY
    % mismatches across all the dimensions.
any(any(any(cfg.vol~=out.vol.data))) % Should only return 1 if there's disagreement somewhere..
    % Where are these mismatches?
idx = find(cfg.vol~=out.vol.data);
fprintf("input: ") % All ones
fprintf("%d, ",cfg.vol(idx))
fprintf("\n")
fprintf("output: ") % Way bigger (?!)
fprintf("%d, ",out.vol.data(idx))
fprintf("\n")
    % And they agree again JUST outside of these entried!
fprintf("They agree just outside:")
cfg.vol(idx(end)+1) 
out.vol.data(idx(end)+1) 