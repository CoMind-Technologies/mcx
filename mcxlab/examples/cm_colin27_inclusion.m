%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tdNIRS simulation in colin27 voxel tissue model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clear cfg;

%% Define simulation parameters
% set seed to make the simulation repeatible
cfg.seed=hex2dec('623F9A9E'); 
% Volume model being used.
VOLUME='brain';
inclusion = 0;
cfg.nphoton=1e4;

load colin27_v3.mat
cfg.vol=colin27;
% Each line below defines the optical parameters for that tissue type:
% mua (absorption) mus (scattering) g (anisotropy) and 'n'
cfg.prop=[  0         0         1.0000    1.0000 % (0) background/air
        0.0190    7.8182    0.8900    1.3700 %(1) scalp
        0.0190    7.8182    0.8900    1.3700 %(2) skull
        0.0040    0.0090    0.8900    1.3700 %(3) csf
        0.0200    9.0000    0.8900    1.3700 %(4) gray matters
        0.0800   40.9000    0.8400    1.3700 %(5) white matters
             0         0    1.0000    1.0000]; %(6) air pockets
         
cfg.srcpos=[75 14 73]; %Source position
cfg.detpos=[75 18 92 5]; %Detector Position
yz_plane = int8(0.5*(cfg.srcpos(1)+cfg.detpos(1))); % Define the planes for plotting
xy_plane = int16(0.5*(cfg.srcpos(3)+cfg.detpos(3)));
% if inclusion
%     incl_origin = [75,35, 85];
%     incl_radius = 2;
%     cfg.shapes = sprintf('{"Shapes":[{"Sphere": {"Tag":7, "O":[%d,%d,%d],"R":%d}}]}',[incl_origin,incl_radius]);
%     cfg.prop = [cfg.prop; 100.0 0 0 1.0000]; %(7) inclusion
%     yz_plane = int8(0.5*(cfg.srcpos(1)+cfg.detpos(1))) % Define the planes for imagine later
%     xy_plane = int8(0.5*(cfg.srcpos(3)+cfg.detpos(3)))
% end

% time-domain simulation parameters
cfg.tstart=0;% Start
cfg.tend=5e-9;% Stop
cfg.tstep=2e-10;% timestep

cfg.srcdir=[0.0 1.0 0.0, 5.0];% define the source DIRECTION
cfg.issrcfrom0 = 0; % Defines the location of the origin: 0 or 1 (0 here).
cfg.isreflect=0; % enable reflection at exterior boundary
cfg.isrefint=1;  % enable reflection at interior boundary too

%% Define OUTPUT settings 
cfg.issavedet=0; % Record partial pathlength of detected photons
cfg.ismomentum=0; % Save photon momentums
cfg.issaveref=1; % Save diffuse reflectance at air voxels.
cfg.issaveexit=1; % Save photon exit postitions & directions.
cfg.replaydet=-1; % Replay photons from all detectors.
cfg.issaveseed=1; % save the initial seed
cfg.maxdetphoton=1e6; % Maximum number of photons to detect.

%% GPU thread configuration
cfg.autopilot=1;% Let MCX decide how many threads are important
cfg.gpuid=1;% 

%% Run Simulation
fprintf('running simulation ...\n');
tic;
%[f2,det2]=mcxlabcl(cfg);
[fluence,detphoton,vol,seed,trajectory] = mcxlabcl(cfg);
toc;

% Everything above here runs fine. The next section throws the "Invalid
% Buffer" error

%% Replay detected photons

% newcfg=cfg;
% newcfg.seed=seed.data;
% newcfg.outputtype='jacobian';
% newcfg.detphotons=detphoton.data;
% newcfg.replaydet=1;
% [flux2, detp2, vol2, seeds2, traj2]=mcxlabcl(newcfg);
% jac=sum(flux2.data,4); 


%% plot the results
figure
 % Plot from LOADED volume data 
% horiz = single(squeeze(colin27(:, :, xy_plane)))';
% saggital = single(squeeze(colin27(yz_plane, :, :)))';
    % Plot from OUTPUT volume data
horiz = single(squeeze(vol.data(:, :, xy_plane)))';
saggital = single(squeeze(vol.data(yz_plane, :, :)))';
    % Plot from INPUT volume data
% horiz = single(squeeze(cfg.vol(:, :, xy_plane)))';
% saggital = single(squeeze(cfg.vol(yz_plane, :, :)))';

subplot(121); % plot x,y plane
set(gca,'YDir','normal') 
% contourf(log10(squeeze(sum(fluence.data(:,:,xy_plane,:),4))'), 12, 'LineStyle','none');
% colorbar;
title('total photon flux (horizontal plane)');
hold on;
imcontour(horiz, 'k');
hold on;
% if inclusion
%     horiz_inclusion = single(vol.data(:, :, incl_origin(3)))';
%     imcontour(horiz_inclusion, 'm');
%     hold on;
% end
% Plot the source & detector
plot(cfg.detpos(1), cfg.detpos(2), 'ko')
hold on;
plot(cfg.srcpos(1), cfg.srcpos(2), 'ro')

subplot(122); % plot y,z plane
set(gca,'YDir','normal') 
% contourf(log10(squeeze(sum(fluence.data(yz_plane,:,:,:),4))'), 12, 'LineStyle','none');
% colorbar;
title('total photon flux (saggital plane)');
hold on;
imcontour(saggital, 'k');
hold on;
% if inclusion
%     saggital_inclusion = single(squeeze(vol.data(incl_origin(1), :, :)))';
%     imcontour(saggital_inclusion, 'm');
%     hold on;
% end

% Plot source and detector
plot(cfg.detpos(2), cfg.detpos(3), 'ko')
hold on;
plot(cfg.srcpos(2), cfg.srcpos(3), 'ro')