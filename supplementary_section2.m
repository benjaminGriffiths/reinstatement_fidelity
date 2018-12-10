%% EEG: Control Analyses
% prepare workspace
clearvars
close all
clc

% define root directory
if ispc;    dir_bids = 'Y:/projects/reinstatement_fidelity/bids_data/';
            dir_tool = 'Y:/projects/general/';
            dir_repos  = 'E:/bjg335/projects/reinstatement_fidelity/bids_data/'; % repository directory
else;       dir_bids = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity';
            dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general';
end

% add subfunctions
addpath([dir_bids,'scripts/subfunctions'])

% define number of subjects
n_subj = 21;

%% Prepare Occipital ROI
% load template grid
load('Y:/projects/general/fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d4mm.mat'); 
template_grid = sourcemodel; clear sourcemodel

% load AAL atlas
mri = ft_read_mri([dir_tool,'fieldtrip-20170319/template/atlas/aal/ROI_MNI_V4.nii']);
 
% reshape mri.anatomy
orig_shape  = size(mri.anatomy);
mri.anatomy = mri.anatomy(:);

% change MRI to binary 'in-occipital' vs. 'out-occipital'
mri.anatomy = mri.anatomy == 5001 | mri.anatomy == 5002 | ... % calcarine L/R 
    mri.anatomy == 5011 | mri.anatomy == 5012 | ... % cuneus L/R 
    mri.anatomy == 5021 | mri.anatomy == 5022 | ... % lingual L/R 
    mri.anatomy == 5101 | mri.anatomy == 5102 | ... % occipital superior L/R 
    mri.anatomy == 5201 | mri.anatomy == 5202 | ... % occipital middle L/R 
    mri.anatomy == 5301 | mri.anatomy == 5302; % occiptial inferior L/R 

% reshape mri.anatomy
mri.anatomy = reshape(mri.anatomy,orig_shape);

% interpolate mri with grid
cfg                 = [];
cfg.parameter       = 'anatomy';
cfg.interpmethod    = 'nearest';
roi                 = ft_sourceinterpolate(cfg,mri,template_grid);

% clean up
clear mri template_grid cfg

%% Get TF of Source Data
% predefine matrices for data
group_erp = nan(n_subj,551);
group_fft = nan(n_subj,275);
fft_freq  = 100*(1:275)/551;
erp_time  = linspace(-1.75,3.75,551);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject handle
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % load
    load([dir_bids,'sourcedata/',subj_handle,'/eeg/',subj_handle,'_task-rf_eeg.mat'])
        
    % limit data to occipital ROI
    for trl = 1 : numel(source.trial)
        source.trial{trl} = source.trial{trl}(find(roi.anatomy(brainGeo.inside==1)),:); %#ok<FNDSB> <- weird bug relating to boolean values!?
    end
    
    % update channel information
    source.label    = source.label(find(roi.anatomy(brainGeo.inside==1))); %#ok<FNDSB>
    
    % get timelocked data
    tml             = ft_timelockanalysis([],source);
    
    % filter
    cfg             = [];
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 15;
    tml             = ft_preprocessing(cfg,tml);
    
    % take absolute of timelocked data (address phase differences between subjs)
    tml.avg = abs(tml.avg);    
    
    % baseline correct timelocked data
    cfg             = [];
    cfg.baseline    = [-0.25 0];
    tml             = ft_timelockbaseline(cfg,tml);
    
    % average over channels
    cfg             = [];
    cfg.avgoverchan = 'yes';
    tml             = ft_selectdata(cfg,tml);
    
    % extract data
    group_erp(subj,:) = tml.avg;
        
    % predefine matrix for subject FFT
    subj_fft = nan([numel(source.trial),size(source.trial{1})]);
    
    % cycle through each trial
    for trl = 1 : numel(source.trial)
        for chan = 1 : numel(source.label)
            subj_fft(trl,chan,:) = abs(fft(source.trial{trl}(chan,:)));
        end
    end
    
    % take first half of output (excluding DC component)
    subj_fft = subj_fft(:,:,2:ceil(size(subj_fft,3)/2));
    
    % add average to group
    group_fft(subj,:) = squeeze(mean(mean(subj_fft,2),1));
    
    % clean workspace
    clear subj_handle source trl cfg tml subj_fft chan
end

% save
save([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-sanity.mat'],'group_erp','group_fft','fft_freq','erp_time')

%% Plot
h=figure('units','centimeters','position',[1 2 17.8 6]); hold on
subplot(1,2,1); hold on
shadedErrorBar(erp_time,mean(abs(group_erp)),sem(abs(group_erp)))
plot([0 0],[0 2000],'k--')
xlabel('Time (s)')
ylabel('Abs. Amp. (µV)'); ylim([0 2000])
xlim([-0.25 1])
title('a','position',[-0.6 2000],'FontSize',14)
set(gca,'tickdir','out','box','off','xtick',-0.25:0.25:1,'ytick',0:400:2000)

subplot(1,2,2); hold on
shadedErrorBar(fft_freq,mean(group_fft),sem(group_fft))
xlabel('Freq. (Hz)')
ylabel('Power (a. u.)'); ylim([0 1200000])
title('b','position',[-10 1200000],'FontSize',14)
set(gca,'tickdir','out','box','off','xtick',0:10:50,'ytick',(0:2:12)*100000)

print(h,[dir_repos(1:end-10),'figures/supplementary2'],'-dtiff','-r300')

%% Run Whole Brain Timelock Analysis
% load template grid
load('Y:/projects/general/fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d4mm.mat'); 
template_grid = sourcemodel; clear sourcemodel

% predefine matrices for data
group_erp = nan(n_subj,23156);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject handle
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % load
    load([dir_bids,'sourcedata/',subj_handle,'/eeg/',subj_handle,'_task-rf_eeg.mat'])
        
    % get timelocked data
    tml             = ft_timelockanalysis([],source);
    
    % filter
    cfg             = [];
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 15;
    tml             = ft_preprocessing(cfg,tml);
    
    % take absolute of timelocked data (address phase differences between subjs)
    tml.avg = abs(tml.avg);    
    
    % baseline correct timelocked data
    cfg             = [];
    cfg.baseline    = [-0.25 0];
    tml             = ft_timelockbaseline(cfg,tml);
    
    % average over channels
    cfg             = [];
    cfg.latency     = [0.1 0.2];
    cfg.avgovertime = 'yes';
    tml             = ft_selectdata(cfg,tml);
    
    % extract data
    group_erp(subj,:) = tml.avg;
        
    % clean workspace
    clear subj_handle source trl cfg tml subj_fft chan
end

% create source data structure
source                  = [];
source.inside           = brainGeo.inside;
source.dim              = brainGeo.dim;
source.pos              = template_grid.pos*10;
source.unit             = 'mm';

% define powspctrm of cluster
source.pow              = nan(size(brainGeo.pos,1),1);     
source.pow(brainGeo.inside==1)	= nanmean(group_erp,1); 

% reshape data to 3D
source.pow              = reshape(source.pow,source.dim);
source.inside           = reshape(source.inside,source.dim);

% add transformation matrix
source.transform        = [1,0,0,-91;
                           0,1,0,-127;
                           0,0,1,-73;
                           0,0,0,1];

% export
cfg = [];
cfg.parameter     = 'pow';               % specify the functional parameter to write to the .nii file
cfg.filename      = [dir_bids,'derivatives/group/eeg/group_task-rf_eeg-erpmap.nii'];  % enter the desired file name
cfg.filetype      = 'nifti';
cfg.coordsys      = 'spm';
cfg.vmpversion    = 2;
cfg.datatype      = 'float';
ft_volumewrite(cfg,source);      % be sure to use your interpolated source data
