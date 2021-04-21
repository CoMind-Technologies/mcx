%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Photonics simulation via MCXLAB in a simple tissue model built up in
% script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clear cfg;
% Below few lines tries to run mcx; failing which, assumes must use mcxcl
try
    gpuinfo=mcxlab('gpuinfo');
catch
    USE_MCXCL=1;
end

%% Define bare minimum simulation configuration (cfg)
cfg.nphoton=1e8; % Number of photons
cfg.vol = uint8(zeros(110,110,110)); % A uniform 'null' zone (0)
cfg.vol(5:105,5:105,5:105)=1; % Material 1 zone
cfg.prop = [0 0 1 1; % Line 1 is for zeros, but these are 'vanishing' points.
            0.0680        1.0800      0.9200       1.3700]; % Line2 for material 1: 
            % absorption, scattering, anisotropy & refractive index.
cfg.srcpos=[80 5 50]; %Source position
cfg.srcdir=[0.0 1.0 0.0, 5.0];% define the source DIRECTION
cfg.tstart=0;% Start time
cfg.tend=8e-9;% Stop time
cfg.tstep=1e-9;% timestep (sec)
%% We can preview our setup:
mcxpreview(cfg)
%% Add somelayers to scatter things
cfg.vol(5:105,10:30,5:105)=2; % Material type 2
cfg.vol(5:105,31:60,5:105)=3; % Material type 3
cfg.prop = [cfg.prop(1:2,:);% Adding properties
    0.0190    7.8182    0.8900    1.3700;%Like Bone
    0.0200    3.8182    0.9000    1.3700]; % highly absorbative
% Include a detector
cfg.detpos = [40,cfg.srcpos(2),cfg.srcpos(3),10]; % [x,y,z,radius]
%% Preview
mcxpreview(cfg)
%% Run a simulation and visualise the data:
cfg.bc='a'; % Impose an absorbative boundary condition
[fluence,detphotons,~,seeds,trajectories]=mcxlab(cfg);
%% Visualise again
horiz = squeeze(cfg.vol(:, :, cfg.srcpos(3)));
figure
subplot(121);
hold on
set(gca,'YDir','normal');
contourf(log10(squeeze(sum(fluence.data(:,:,cfg.srcpos(3),:),4))), 12, 'LineStyle','none');
hold on
imcontour(horiz,'k')
hold on
plot([cfg.detpos(2),cfg.detpos(2)],[cfg.detpos(1)-10,cfg.detpos(1)+10], '-r','LineWidth',2);
title('total photon flux (horizontal plane)');
xlabel('Y')
ylabel('X')
colorbar;

subplot(122)
histogram(detphotons.data(2,:))
set(gca,'Yscale','log');
xlabel("Momentum at detection [a.u.]")
ylabel("Count")
%% Now investigate replay
cfg.seed=seeds.data; % Rertrieve the photon seeds
cfg.outputtype='jacobian'; % Request the Jacobian
cfg.detphotons=detphotons.data; % We want the same detected photons
cfg.replaydet=1;
[newFluence, newDetphotons,~,~,newTrajectories]=mcxlab(cfg);
jacobian=sum(newFluence.data,4); % The jacobian, summed over time gates
%% And visualise
figure(1)
subplot(121)
hold on
contourf(log10(squeeze(sum(newFluence.data(:,:,cfg.srcpos(3),:),4))));
hold on
imcontour(horiz,'k')
title('Detected photon fluence in source (y,x)-plane');
xlabel("Y")
ylabel("X")
colorbar

subplot(122)
hold on;
xt=1e9*((cfg.tstart:cfg.tstep:cfg.tend-cfg.tstep)+0.5*cfg.tstep);
pos = cfg.srcpos(2);
p_end=40;
subset=[];
for i=pos:pos+p_end
    plt=semilogy(xt,squeeze(newFluence.data(cfg.detpos(1),i,cfg.detpos(3),:)),'color',[.8 i/(pos+1.1*p_end) i/(pos+p_end)],'LineWidth',2);
    if (i==pos)|(i==pos+p_end)
        subset=[subset,plt];
    end
end
legend(subset,{'shallow','deep'})
xlabel('time (ns)')
ylabel('fluence (1/mm^2/s)')
set(gca,'yscale','log');
title('Detected photon time point spread functions (TPSF)');
