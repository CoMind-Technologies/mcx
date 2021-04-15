%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A"how it works" script, making evident some functionality of MCX
% Idea is to show a few things:
%   1) How the detector location effects the "pre-processing" of the volume
%   shape.
%   2) how adding an inclusion effects pre-processing of volume shape.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clear cfg;
try % Allows for use of MXCLABCL is on a non-NVIDIA machine.
    gpuinfo=mcxlab('gpuinfo');
catch
    USE_MCXCL=1;
end
%% Upload a default volume
load colin27_v3.mat
cfg.vol=colin27; % The default volume is 
% Each line below defines the optical parameters for that tissue type:
% mua (absorption) mus (scattering) g (anisotropy) and 'n'
cfg.prop=[  0         0         1.0000    1.0000 % (0) background/air
        0.0190    7.8182    0.8900    1.3700 %(1) scalp
        0.0190    7.8182    0.8900    1.3700 %(2) skull
        0.0040    0.0090    0.8900    1.3700 %(3) csf
        0.0200    9.0000    0.8900    1.3700 %(4) gray matters
        0.0800   40.9000    0.8400    1.3700 %(5) white matters
             0         0    1.0000    1.0000]; %(6) air pockets
% %         
cfg.srcpos=[75 14 73]; %Source position
cfg.srcdir=[0.0 1.0 0.0, 5.0];% define the source DIRECTION
% time-domain simulation parameters
cfg.tstart=0;% Start
cfg.tend=5e-9;% Stop
cfg.tstep=2e-10;% timestep
% Above are all essential. Next is just a choice.
cfg.issrcfrom0 = 0; % Defines the location of the origin: 0 or 1 (0 here).
cfg.nphoton=0; % Number of photons
%% Visualise it:
mcxpreview(cfg)
%% Now "process" it through mcx:
[fluence,detphoton,vol] = mcxlab(cfg);
%% And visualise output...
out=cfg;
out.vol=vol.data;
mcxpreview(out);
% Above should work, and look the same!
%% Now include an inclusion
incl_origin = [75 25 92];
incl_radius = 2;
cfg.shapes = [sprintf('{"Shapes":[{"Sphere": {"Tag":7, "O":[%d,%d,%d],"R":%d}}]}',[incl_origin,incl_radius])];
cfg.prop = [cfg.prop(1:7,:); 0.0800   40.9000    0.8400    1.3700]; %(7) inclusion
%yz_plane = int8(0.5*(cfg.srcpos(1)+cfg.detpos(1))) % Define the planes for imagine later
%xy_plane = int8(0.5*(cfg.srcpos(3)+cfg.detpos(3)))
% And visualise:
mcxpreview(cfg);
%% Now "process" it through mcx:
[fluence,detphoton,vol] = mcxlab(cfg);
%% And visualise output...
out.vol=vol.data;
mcxpreview(out);
% Now what's weird here is ONLY thr inclusion remains?!
%% Now include a detector
cfg.detpos=[75 18 92 5]; %Detector Position
yz_plane = int8(0.5*(cfg.srcpos(1)+cfg.detpos(1))); % Define the planes for plotting
xy_plane = int16(0.5*(cfg.srcpos(3)+cfg.detpos(3)));

% cfg.isreflect=0; % enable reflection at exterior boundary
% cfg.isrefint=1;  % enable reflection at interior boundary too

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
[fluence,detphoton,vol,seed,trajectory] = mcxlab(cfg);
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