%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Builds a 'default' slab model and a set of identical volumes with
% inclusions.
% Saves a default configuration file and (separately) the inclusion volumes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clear cfg;
try
    gpuinfo=mcxlab('gpuinfo');
catch
    USE_MCXCL=1;
end
%% Setup
% All numbers in mm. Will be converted to volume via cfg.unitmm
% (effective spatial resolution)
cfg.unitinmm=0.5; % voxels are (0.5 mm)^3 volumes
dims = [50,50,100]; % x,y,z dimensions in mm
cfg.vol=zeros(dims(1)/cfg.unitinmm,dims(2)/cfg.unitinmm,dims(3)/cfg.unitinmm,'uint8');
cfg.dims=size(cfg.vol); % Could cause a problem?
% 
surface_z = double((dims(3)-30)/cfg.unitinmm+1);
cfg.vol(:,:,surface_z:end)=5; % 3cm of air at surface
cfg.vol(:,:,end-35/cfg.unitinmm+1:end-30/cfg.unitinmm)=1; % 5 mm of scalp
cfg.vol(:,:,end-42/cfg.unitinmm+1:end-35/cfg.unitinmm)=2; % 7 mm of skull
cfg.vol(:,:,end-47/cfg.unitinmm+1:end-42/cfg.unitinmm)=3; % 1 mm of csf
cfg.vol(:,:,1:end-47/cfg.unitinmm)=4; % etc. is 'grey-matter'

%          % absorp | red. scat | anisotropy | rarefrac. index 
cfg.prop=[  0         0         1.0000    1.0000 % (0) background/air
            0.0190    7.8182    0.8900    1.3700 % (1) scalp
            0.0190    7.8182    0.8900    1.3700 % (2) skull
            0.0040    0.0090    0.8900    1.3700 % (3) csf
            0.0200    9.0000    0.8900    1.3700 % (4) gray matters
%            0.0800   40.9000    0.8400    1.3700 % (-) white matters
                 0         0    1.0000    1.0000]; % (5) air pockets
cfg.srcpos=[0.5*cfg.dims(1) 0.5*cfg.dims(2) surface_z]; % Source position
cfg.srcdir=[0 0 -1]; % Source direction
% Detector details
det_distance=30; % distance from source, in mm
det_radius=3; % detector radius, in mm
cfg.detpos=[cfg.srcpos(1)-sqrt(det_distance/cfg.unitinmm) cfg.srcpos(2)-sqrt(det_distance/cfg.unitinmm) cfg.srcpos(3) det_radius/cfg.unitinmm % Detector 1
            cfg.srcpos(1)+sqrt(det_distance/cfg.unitinmm) cfg.srcpos(2)-sqrt(det_distance/cfg.unitinmm) cfg.srcpos(3) det_radius/cfg.unitinmm % Detector 2
            cfg.srcpos(1)-sqrt(det_distance/cfg.unitinmm) cfg.srcpos(2)+sqrt(det_distance/cfg.unitinmm) cfg.srcpos(3) det_radius/cfg.unitinmm % Detector 2
            cfg.srcpos(1)+sqrt(det_distance/cfg.unitinmm) cfg.srcpos(2)+sqrt(det_distance/cfg.unitinmm) cfg.srcpos(3) det_radius/cfg.unitinmm]; % Detector 3
%mcxpreview(cfg)
save('../configs/fourLayer_slab.mat','cfg')
%% Creating set of inclusions
% inclusion ranges differ by depth and size
rng('default') % Default the rng
inclusion_depths =[14 18 22 26 30 34 38 42]; % in mm, from surface
inclusion_widths = [0.5 1 2 4];% in mm.
% create inclusion indices, in voxel space
inclusion_set=[];
for d = inclusion_depths
    for w = inclusion_widths
        x0 = cfg.srcpos(1)-floor(0.5*w/cfg.unitinmm);
        x1 = cfg.srcpos(1)+ceil(0.5*w/cfg.unitinmm);
        y0 = cfg.srcpos(2)-floor(0.5*w/cfg.unitinmm);
        y1 = cfg.srcpos(2)+ceil(0.5*w/cfg.unitinmm);
        z1 = surface_z-(d/cfg.unitinmm);
        z0 = z1-(w/cfg.unitinmm);
        incl.xyz=meshgrid(x0:x1,y0:y1,z0:z1);
        incl.idx=sub2ind(cfg.dims,incl.xyz);
        inclusion_set=[inclusion_set incl];
    end
end
% Properties of inlcusion (haemoglobin?)
incl_prop = [0.0800    9.0000    0.8900    1.3700];
%% Ammend with simulation properties
cfg.nphoton=1e7;
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