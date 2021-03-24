%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tdNIRS simulation in colin27 voxel tissue model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clear cfg;

%% preparing the input data
% set seed to make the simulation repeatible
cfg.seed=hex2dec('623F9A9E'); 

VOLUME='brain';

cfg.nphoton=1e8;

if strcmp(VOLUME, 'brain')
    load colin27_v3.mat
    cfg.vol=colin27;
    cfg.prop=[  0         0         1.0000    1.0000 % background/air
            0.0190    7.8182    0.8900    1.3700 % scalp
            0.0190    7.8182    0.8900    1.3700 % skull
            0.0040    0.0090    0.8900    1.3700 % csf
            0.0200    9.0000    0.8900    1.3700 % gray matters
            0.0800   40.9000    0.8400    1.3700 % white matters
                 0         0    1.0000    1.0000]; % air pockets
    cfg.srcpos=[75 14 73];
    cfg.detpos=[75 18 92 5];
    
elseif strcmp(VOLUME, 'cube')
    dim=60;
    [xi,yi,zi]=meshgrid(1:dim,1:dim,1:dim);
    dist=(xi-30).^2+(yi-30).^2+(zi-30).^2;
    cfg.vol=ones(size(xi));
    cfg.vol(dist<100)=2;
    cfg.vol=uint8(cfg.vol);
    cfg.srcpos=[30,30,0]+1;
    cfg.prop=[0 0 1 1          % medium 0: the environment
       0.002 1.0 0.01 1.37     % medium 1: cube
       0.050 0.5 0.01 1.37];   % medium 2: spherical inclusion
end

% define the source position

cfg.srcdir=[0.0 1.0 0.0, 5.0];
%cfg.srcdir=[0.0 0.0 -1.0 5.0];
%cfg.srctype='gaussian';
%cfg.srcparam1=[10 0 0 0];
%cfg.srcparam2=[0 0 0 0];

% 1 = first voxel is [0 0 0]
cfg.issrcfrom0 = 0;

% time-domain simulation parameters
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=2e-10;

% GPU thread configuration
cfg.autopilot=1;
cfg.gpuid=1;

%% running simulation with boundary reflection enabled
fprintf('running simulation ...\n');
cfg.isreflect=1; % enable reflection at exterior boundary
cfg.isrefint=1;  % enable reflection at interior boundary too
cfg.issavedet=1; % enable recording partial pathlength of detected photons
cfg.ismomentum=1
cfg.issaveref=1
cfg.issaveexit=1;
cfg.replaydet=1

tic;
%[f2,det2]=mcxlab(cfg);
[fluence,detphoton,vol,seed,trajectory] = mcxlab(cfg);
toc;

%%
% Replay detected photons

newcfg=cfg;
newcfg.seed=seed.data;
newcfg.outputtype='jacobian';
newcfg.detphotons=detphoton.data;
newcfg.replaydet=1
[flux2, detp2, vol2, seeds2, traj2]=mcxlab(newcfg);
jac=sum(flux2.data,4); 


%% plot the results
figure

subplot(231); % plot x,y plane
horiz = single(colin27(:, :, cfg.detpos(3)))';
set(gca,'YDir','normal') 
contourf(log10(squeeze(sum(fluence.data(:,:,cfg.detpos(3),:),4))'), 12, 'LineStyle','none');
title('total photon flux (horizontal plane)');
colorbar;
hold on;
imcontour(horiz, 'k');
hold on;
plot(cfg.detpos(1), cfg.detpos(2), 'ko')
hold on;
plot(cfg.srcpos(1), cfg.srcpos(2), 'ro')

subplot(232); % plot y,z plane
hold on;
saggital = single(squeeze(colin27(cfg.detpos(1), :, :)))';
set(gca,'YDir','normal') 
contourf(log10(squeeze(sum(fluence.data(cfg.detpos(1),:,:,:),4))'), 12, 'LineStyle','none');
colorbar;
title('total photon flux (saggital plane)');
hold on;
imcontour(saggital, 'k');
hold on;
plot(cfg.detpos(2), cfg.detpos(3), 'ko')
hold on;
plot(cfg.srcpos(2), cfg.srcpos(3), 'ro')

subplot(233);
hold on;
xi=1e9*((cfg.tstart:cfg.tstep:cfg.tend-cfg.tstep)+0.5*cfg.tstep);
pos = cfg.srcpos(2);
for i=pos:pos+10
    semilogy(xi,squeeze(flux2.data(cfg.detpos(1),i,cfg.detpos(3),:)),'color',[1-(i-pos)/pos 1-(i-pos)/pos 1]);
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
hold on;
plot(cfg.detpos(2), cfg.detpos(3), 'ko')
hold on;
plot(cfg.srcpos(2), cfg.srcpos(3), 'ro')
colorbar;

subplot(236);
hold on;
[hs1,c1]=hist(detphoton.ppath(find(detphoton.ppath(:,1)),1),200);
[hs2,c2]=hist(detphoton.ppath(find(detphoton.ppath(:,2)),2),200);
[hs3,c3]=hist(detphoton.ppath(find(detphoton.ppath(:,3)),3),200);
[hs4,c4]=hist(detphoton.ppath(find(detphoton.ppath(:,4)),4),200);
[hs5,c5]=hist(detphoton.ppath(find(detphoton.ppath(:,5)),5),200);
[hs6,c6]=hist(detphoton.ppath(find(detphoton.ppath(:,6)),6),200);
bar(c1,hs1,'edgecolor','none','facecolor','r');
bar(c2,hs2,'edgecolor','none','facecolor','y');
bar(c3,hs3,'edgecolor','none','facecolor','k');
bar(c4,hs4,'edgecolor','none','facecolor','b');
h=bar(c5,hs5,'edgecolor','none','facecolor','k');
h.FaceAlpha = 0.2;
bar(c6,hs6,'edgecolor','none','facecolor','c');
legend('air','scalp','skull','csf','gray','white');
xlabel('partial pathlength (mm)');
title(sprintf('detected %d photons',size(detphoton.ppath,1)));
