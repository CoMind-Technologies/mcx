%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tdNIRS simulations in head-tissue slab model with varying inclusions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clear cfg;
try
    gpuinfo=mcxlab('gpuinfo');
catch
    USE_MCXCL=1;
end
%% Load prepared config file
load colin27_mc_slab.mat;
%% Ammend with simulation properties
cfg.nphoton=1e8;
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=2e-10;

fprintf('running simulation ...\n');
cfg.isreflect=0; % enable reflection at exterior boundary
cfg.isrefint=1;  % enable reflection at interior boundary too
cfg.issavedet=1; % enable recording partial pathlength of detected photons
cfg.ismomentum=1;
cfg.issaveref=1;
cfg.issaveexit=1;
% GPU Config
cfg.autopilot=1;
cfg.gpuid=1;
%% Running simulation
tic;
[fluence,detphoton,vol,seed,trajectory] = mcxlab(cfg);
toc;
%% Get Jacobian
% cfg.detphoton=detphoton.data;
cfg.seed=seed.data;
cfg.detphotons = detphoton.data;
cfg.outputtype='jacobian';
cfg.replaydet=0;
%% 
[flux2, detp2, vol2, seeds2,traj2]=mcxlab(cfg);
%% plot the results 1652
figure
subplot(231); % plot x,y plane
horiz = squeeze(vol.data(:, :, cfg.detpos(3)))';
horiz(horiz>1000) = 1;
horiz = uint8(horiz);
set(gca,'YDir','normal');
contourf(log10(squeeze(sum(fluence.data(:,:,cfg.detpos(3),:),4))'), 12, 'LineStyle','none');
title('total photon flux (horizontal plane)');
xlabel('X')
ylabel('Y')
colorbar;
hold on;
imcontour(horiz, 'k');
hold on;
plot(cfg.detpos(1), cfg.detpos(2), 'ko');
hold on;
plot(cfg.srcpos(1), cfg.srcpos(2), 'ro');

subplot(232); % plot x,z plane
hold on;
vertical = squeeze(vol.data(:,cfg.detpos(2), :))';
vertical(vertical>1000) = 1;
vertical= uint8(vertical);
set(gca,'YDir','normal') 
% Plot the flux through the vertical plane half-way between S & D
contourf(log10(squeeze(sum(fluence.data(:,0.5*(cfg.detpos(2)+cfg.srcpos(2)),:,:),4))'), 12, 'LineStyle','none');
colorbar;
title('total photon flux (vertical plane)');
xlabel('X')
ylabel('Z')
hold on;
imcontour(vertical, 'k');
hold on;
plot(cfg.detpos(1), cfg.detpos(3), 'ko')
hold on;
plot(cfg.srcpos(1), cfg.srcpos(3), 'ro')
colorbar;

subplot(233); % TPSF
hold on;
xt=1e9*((cfg.tstart:cfg.tstep:cfg.tend-cfg.tstep)+0.5*cfg.tstep);
pos = cfg.srcpos(2);
for i=pos:pos+10
    semilogy(xt,squeeze(flux2.data(cfg.detpos(1),i,cfg.detpos(3),:)),'color',[1-(i-pos)/pos 1-(i-pos)/pos 1]);
end
xlabel('time (ns)')
ylabel('flux (1/mm^2/s)')
set(gca,'yscale','log');
title('detected photon time point spread functions (TPSF)');

subplot(234); % plot x,y plane
hold on;
contourf(log10(squeeze(sum(flux2.data(:,:,cfg.detpos(3),:),4))'), 12, 'LineStyle','none');
title('detected photon flux (horizontal plane)');
hold on;
imcontour(horiz, 'k');
xlabel('X')
ylabel('Y')
hold on;
plot(cfg.detpos(1), cfg.detpos(2), 'ko')
hold on;
plot(cfg.srcpos(1), cfg.srcpos(2), 'ro')
colorbar;

subplot(235); % plot x,y plane
hold on;
contourf(log10(squeeze(sum(flux2.data(:,0.5*(cfg.detpos(2)+cfg.srcpos(2)),:,:),4))'), 12, 'LineStyle','none');
title('detected photon flux (vertical plane)');
hold on;
imcontour(vertical, 'k');
xlabel('Y')
ylabel('Z')
colorbar;
hold on;
plot(cfg.detpos(1), cfg.detpos(3), 'ko')
hold on;
plot(cfg.srcpos(1), cfg.srcpos(3), 'ro')
colorbar;

subplot(236); % plot momenta
hold on;
[hs1,c1]=hist(detphoton.ppath(find(detphoton.ppath(:,1)),1),200);
bar(c1,hs1,'edgecolor','none','facecolor','r');
xlabel('partial pathlength (mm)');
title(sprintf('detected %d photons',size(detphoton.ppath,1)));