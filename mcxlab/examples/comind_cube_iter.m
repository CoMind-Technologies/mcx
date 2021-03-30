%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tdNIRS simulation in colin27 voxel tissue model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clear cfg;

%% preparing the input data
% set seed to make the simulation repeatible
cfg.seed=hex2dec('623F9A9E'); 

cfg.nphoton=1e7;


dim=60;
[xi,yi,zi]=meshgrid(1:dim,1:dim,1:dim);
%dist=(xi-30).^2+(yi-30).^2+(zi-30).^2;
cfg.vol=ones(size(xi));
%cfg.vol(dist<100)=2;
cfg.vol=uint8(cfg.vol);
cfg.srcpos=[0,30,15];
cfg.detpos=[0,30,45,5];

%[mua, mus, g, n]
cfg.prop=[0 0 1 1          % medium 0: the environment
   0.0680    1.0800    0.9200    1.3700];     % medium 1: gray matter
   %0.08 0.5 0.01 1.37];   % medium 2: spherical inclusion
   
   
src_dir = [0.0 1.0 0.0];
cfg.srcdir = src_dir/norm(src_dir);
cfg.src_dir(4) = 5.0;
cfg.srctype='gaussian';
cfg.srcparam1 = [5 0 0 0];
cfg.srcparam2 = [0 0 0 0];

% 1 = first voxel is [0 0 0]
cfg.issrcfrom0 = 0;

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

all_dat = zeros(size(src_dirs, 1),1);
for n=1:size(src_dirs, 1)
    
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
    
    all_dat(n) = size(detphoton.detid, 1);
end

%%

c=zeros(5,1);
for n=1:5; 
    c(n)=rad2deg(atan(src_dirs(n,3)/src_dirs(n,2)));
end;

plot(c, all_dat);
xlabel('launch angle');
ylabel('number detected photons');