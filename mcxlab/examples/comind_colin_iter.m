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

load colin27_v3.mat;
cfg.vol=colin27;
% 840nm [mua, mus, g, n]
cfg.prop=[0         0         1.0000    1.0000 % background/air
        0.0604    2.5205    0.9200    1.3700 % scalp
        0.0922    4.7953    0.9200    1.3700 % skull
        0.0040    0.0090    0.9200    1.3700 % csf
        0.0680    1.0800    0.9200    1.3700 % gray matter
        0.0681    1.1894    0.9200    1.3700 % white matter
             0         0    1.0000    1.0000]; % air pockets
cfg.srcpos=[75 8 68];
cfg.detpos=[75 14 92 3];

% define the source position

%src_dir = [0.0 1.0 0.0];
%cfg.srcdir = src_dir/norm(src_dir);
%cfg.src_dir(4) = 5.0;
cfg.srctype = 'gaussian';
cfg.srcparam1 = [5 0 0 0];
cfg.srcparam2 = [0 0 0 0];

% 1 = first voxel is [0 0 0]
cfg.issrcfrom0 = 1;

% time-domain simulation parameters
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=2e-10;

% GPU thread configuration
cfg.autopilot=1;
cfg.gpuid=1;

fprintf('running simulation ...\n');
cfg.isreflect=0; % enable reflection at exterior boundary
cfg.isrefint=1;  % enable reflection at interior boundary too
cfg.issavedet=1; % enable recording partial pathlength of detected photons
cfg.ismomentum=1;
cfg.issaveref=1;
cfg.issaveexit=1;
cfg.replaydet=1;

src_dirs = [[0.0 0.0 1.0]; [0.0 0.5 0.5]; [0.0 1.0 0.0]; 
    [0.0 0.5 -0.5]; [0.0 0.0 -1.0]];

all_dat = zeros(size(src_dirs, 1),2);
for n=1:size(src_dirs, 1)
    
    %cfg.detpos=[75 14 92 det_radii(1,n)];
    src_dir = src_dirs(n,:);
    cfg.srcdir = src_dir/norm(src_dir);
    cfg.src_dir(4) = 5.0;

    tic;
    [fluence,detphoton,vol,seed,trajectory] = mcxlab(cfg);
    toc;

    newcfg=cfg;
    newcfg.seed=seed.data;
    newcfg.outputtype='jacobian';
    newcfg.detphotons=detphoton.data;
    newcfg.replaydet=1;
    [flux2, detp2, vol2, seeds2, traj2]=mcxlab(newcfg);
    jac=sum(flux2.data,4); 

    gray_photons = detphoton.ppath(find(detphoton.ppath(:,4)),4);
    gray_total_path = sum(gray_photons);
    all_total_path = sum(sum(detphoton.ppath, 1));
    perc_gray = size(gray_photons,1) / size(detphoton.detid, 1);
    all_dat(n,1) = perc_gray;
    all_dat(n,2) = gray_total_path / all_total_path;
end

%%
c=zeros(5,1);
for n=1:5; 
    c(n)=rad2deg(atan(-src_dirs(n,3)/src_dirs(n,2)));
end;

%all_dat_norm(:,2) = all_dat(:,2) / max(all_dat(:,2));

plot(c, all_dat);
xlabel('launch angle');
ylabel('fraction photons entering grey matter');
legend('fraction photons entering grey matter', 'fraction of total pathlength in gray matter');