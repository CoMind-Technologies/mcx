clear all;
clear cfg;
try
    gpuinfo=mcxlab('gpuinfo');
catch
    USE_MCXCL=1;
end
%% load volume and inclusions
load colin27_v3.mat
load colin27_csf_500vol.mat

%% Create a set of config files, one for each inclusion
