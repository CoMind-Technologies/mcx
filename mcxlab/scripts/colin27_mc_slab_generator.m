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
% Uncomment the below line to preview the tissue
% mcxpreview(cfg)
%% save the config
save('../configs/colin27_mc_slab.mat','cfg')