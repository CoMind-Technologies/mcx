.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% The above is telling us that  all 'skull' tissue is in voxels labelled '2', so they can be identified...
dims = size(cfg.vol);
scalp.idx = find(cfg.vol==1);
[scalp.x,scalp.y,scalp.z] = ind2sub(dims,scalp.idx);

% Pick a point in the tissue type you like. This can be done by pulling a 'row' of coords
% Done here for cerebral spinal fluid, making 100 inclusions:
s = RandStream('mlfg6331_64');
% Randomly pick a set of volume indices
vox_volume=500;
csf.inclusions=[];
i=0;
while (length(csf.inclusions)<100) && i<10000
    idx = randi(s,length(csf.idx));
    proposal = make_inclusion([csf.x(idx),csf.y(idx),csf.z(idx)],vox_volume,cfg.vol);
    if length(proposal.indices)==vox_volume
        csf.inclusions =[csf.inclusions, proposal];
        fprintf("Inclusion added. Total: %d \n",length(csf.inclusions));
    end
   i=i+1; 
end
% This is essentially all we need; an inclusion needs to add its optical properties  to the
% "cfg.prop" object, but we need not make a whole new VOLUME object for
% each inclusion.
%% Save generated inclusion files
csf_inclusions = csf.inclusions
save(sprintf('./inclusions/csf_%dvol.mat',vox_volume),'csf_inclusions');
clear csf_inclusions;
%% To visualise, we require some bare-minima
cfg.srcpos=[75 14 73]; %Source position
cfg.srcdir=[0.0 1.0 0.0, 5.0];% define the source DIRECTION
% time-domain simulation parameters
cfg.tstart=0;% Start
cfg.tend=5e-9;% Stop
cfg.tstep=cfg.tend;% timestep can be same as duration for no photons.
% Above are all essential. Next is just a choice.
cfg.issrcfrom0 = 0; % Defines the location of the origin: 0 or 1 (0 here).
cfg.nphoton=0; % Number of photons can be zero for generating.
% Each line below defines the optical parameters for that tissue type:
% mua (absorption) mus (scattering) g (anisotropy) and 'n'
cfg.prop=[  0         0         1.0000    1.0000 % (0) background/air
        0.0190    7.8182    0.8900    1.3700 %(1) scalp
        0.0190    7.8182    0.8900    1.3700 %(2) skull
        0.0040    0.0090    0.8900    1.3700 %(3) csf
        0.0200    9.0000    0.8900    1.3700 %(4) gray matters
        0.0800   40.9000    0.8400    1.3700 %(5) white matters
             0         0    1.0000    1.0000 %(6) air pockets
          0.11         1.   0.6       1.37]; %(7) inclusions (blood?)
%%  Try find the larget ones?
maxes= [0 0 0 0];
idx = [0 0 0 0];
for i=1:length(csf.inclusions)
    max_z = max(csf.inclusions(i).voxels(:,3));
    if maxes(end)<max_z
        maxes(1:3)=maxes(2:4);
        idx(1:3)=idx(2:4);
        maxes(end) = max_z;
        idx(end)=i;
    end
end
%% First Row
figure(1)
hold on
cfg.vol(csf.inclusions(idx(3)).indices)=7;
title(sprintf('Inclusion near (%d,%d,%d)',csf.inclusions(idx(3)).voxels(1,:)));
hold on;
mcxpreview(cfg);
mcxplotvol(cfg.vol);
%% Second Row
cfg.vol=colin27;
cfg.vol(csf.inclusions(idx(4)).indices)=7;
figure(2)
hold on
title(sprintf('Inclusion near (%d,%d,%d)',csf.inclusions(idx(4)).voxels(1,:)));
hold on;
mcxpreview(cfg);
mcxplotvol(cfg.vol);
%% FUNCTION DEFINITIONS:
function inclusion = make_inclusion(pi,n_voxels,volume)
% This function starts at a point, a voxel, and performs a random walk
% through the volume, starting at 
% recording all voxels visited that satisfy 2 conditions:
% pi is the initial point (a vector of 2 usigned integers)
dims = size(volume);
inclusion.voxels =[pi];
loc=pi;
tissue_type = volume(pi(1),pi(2),pi(3));
s=0;
while (length(inclusion.voxels)<n_voxels) && (s<100*n_voxels) % only try for 100 times the volume size.
    proposals = inclusion.voxels+randi(2,3,1)'-1; % Propose a new voxel neighbour to EVERY member.
    new_rows = ~ismember(proposals,inclusion.voxels,'rows'); % see who's a member 
    proposals=proposals(new_rows,:); % Ditch anything that is alread in the set.
    in_tissue=(volume(sub2ind(dims,proposals(:,1),proposals(:,2),proposals(:,3)))==tissue_type);
    inclusion.voxels=[inclusion.voxels; proposals(in_tissue,:)];
    s=s+1;
end
if length(inclusion.voxels)>n_voxels
    inclusion.voxels=inclusion.voxels(1:n_voxels,:);
end
inclusion.indices = sub2ind(size(volume),inclusion.voxels(:,1),inclusion.voxels(:,2),inclusion.voxels(:,3));
end