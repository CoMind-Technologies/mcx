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
%% Preparing the volume
% Load colin 27
load colin27_v3.mat;
% Remove a 'cube' above right ear
% This core volume is simply 60mm (width)x80cm (depth) x 25cm (height)
cfg.vol= colin27(120:180,70:130,65:90);
cfg.prop=[  0         0         1.0000    1.0000 % (0) background/air
        0.0190    7.8182    0.8900    1.3700 %(1) scalp
        0.0190    7.8182    0.8900    1.3700 %(2) skull
        0.0040    0.0090    0.8900    1.3700 %(3) csf
        0.0200    9.0000    0.8900    1.3700 %(4) gray matters
        0.0800   40.9000    0.8400    1.3700 %(5) white matters
             0         0    1.0000    1.0000]; %(6) air pockets
% Create 
dims = size(cfg.vol);
cfg.srcpos=[0,uint8(0.2*dims(2)),uint8(0.5*dims(3))];
% To place the source, we need to find the outmost (largest x) non-zero cell
cfg.srcpos(1)=find(cfg.vol(:,cfg.srcpos(2),cfg.srcpos(3)),1,'last'); 
cfg.srcdir = [-1.0 0. 0.0];
% Create a detector
cfg.detpos=[0,ceil(0.75*dims(2)),double(cfg.srcpos(3)),2];
cfg.detpos(1)=find(cfg.vol(:,cfg.detpos(2),cfg.detpos(3)),1,'last'); 
%% Preview
mcxpreview(cfg)
%% Simulation parameters
% time-domain simulation parameters
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
cfg.replaydet=1;
% GPU Config
cfg.autopilot=1;
cfg.gpuid=1;
%% Running simulation
tic;
[fluence,detphoton,~,seed,trajectory] = mcxlab(cfg);
toc;
%% Setting up sets of cfg's
rng('default') % Default the rng
numRuns = 100;
inclusion_xyz = meshgrid(18:22,28:32,12:16);
incl_prop = [0.11         1.   0.6       1.37];
cfg_set=[];
for i = 1:3
    cfg_ = cfg;% New spec.
    cfg_.vol(inclusion_xyz)=7;
    cfg_.prop=[cfg.prop; incl_prop];
    cfg_.prop(end,1)=(0.5+rand(1))*incl_prop(1);
    cfg_.replaydet=1;
    cfg_.seed=seed.data;
    cfg_.outputtype='jacobian';
    cfg_.detphotons=detphoton.data;
    cfg_set=[cfg_set,cfg_];
end
%% Run big sim
[flux2, detp2, vol2, seeds2,traj2]=mcxlab(cfg_set);
%% plot the results
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

subplot(232); % plot y,z plane
hold on;
saggital = squeeze(vol.data(cfg.detpos(1), :, :))';
saggital(saggital>1000) = 1;
saggital = uint8(saggital);
set(gca,'YDir','normal') 
contourf(log10(squeeze(sum(fluence.data(cfg.detpos(1),:,:,:),4))'), 12, 'LineStyle','none');
colorbar;
title('total photon flux (saggital plane)');
xlabel('Y')
ylabel('Z')
hold on;
imcontour(saggital, 'k');
hold on;
plot(cfg.detpos(2), cfg.detpos(3), 'ko')
hold on;
plot(cfg.srcpos(2), cfg.srcpos(3), 'ro')

subplot(233);
hold on;
xt=1e9*((cfg.tstart:cfg.tstep:cfg.tend-cfg.tstep)+0.5*cfg.tstep);
pos = cfg.srcpos(2);
for i=pos:pos+10
    semilogy(xt,squeeze(flux2.data(cfg.detpos(1),i,cfg.detpos(3),:)),'color',[1-(i-pos)/pos 1-(i-pos)/pos 1]);
end
xlabel('time (ns)')
ylabel('flux (1/mm^2/s)')
set(gca,'yscale','log');
title('detect photon time point spread functions (TPSF)');

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

subplot(235); % plot y,z plane
hold on;
contourf(log10(squeeze(sum(flux2.data(cfg.detpos(1),:,:,:),4))'), 12, 'LineStyle','none');
title('detected photon flux (saggital plane)');
hold on;
imcontour(saggital, 'k');
xlabel('Y')
ylabel('Z')
hold on;
plot(cfg.detpos(2), cfg.detpos(3), 'ko')
hold on;
plot(cfg.srcpos(2), cfg.srcpos(3), 'ro')
colorbar;

subplot(236);
hold on;
[hs1,c1]=hist(detphoton.ppath(find(detphoton.ppath(:,1)),1),200);
%[hs2,c2]=hist(detphoton.ppath(find(detphoton.ppath(:,2)),2),200);
bar(c1,hs1,'edgecolor','none','facecolor','r');
%bar(c2,hs2,'edgecolor','none','facecolor','k');
h.FaceAlpha = 0.2;
%legend('volume', 'inclusion');
xlabel('partial pathlength (mm)');
title(sprintf('detected %d photons',size(detphoton.ppath,1)));