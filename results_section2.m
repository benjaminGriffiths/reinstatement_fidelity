%% EEG: Source Analysis
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

%% Get TF of Source Data
% cycle through each subject
for subj = 1 : n_subj
    
    % define subject handle
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_bids,'sourcedata/',subj_handle,'/eeg/'];
    
    % load
    load([dir_bids,'sourcedata/',subj_handle,'/eeg/',subj_handle,'_task-rf_eeg.mat'])
    
    % get time-frequency representation of data
    cfg             = [];
    cfg.keeptrials  = 'yes';
    cfg.method      = 'wavelet';
    cfg.width       = 6;
    cfg.output      = 'pow';
    cfg.toi         = -1:0.1:3;
    cfg.foi         = 8:2:30; % 100hz sampling
    cfg.pad         = 'nextpow2';
    freq            = ft_freqanalysis(cfg, source); clear source

    % convert powspctrm to single
    freq.powspctrm = single(freq.powspctrm);
    
    % z-transform
    avg_pow = repmat(nanmean(nanmean(freq.powspctrm,4),1),[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
    std_pow = repmat(nanstd(nanmean(freq.powspctrm,4),[],1),[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
    freq.powspctrm = (freq.powspctrm - avg_pow) ./ std_pow; clear avg_pow std_pow
    
    % average over time and frequency
    cfg             = [];
    cfg.avgovertime = 'yes';
    cfg.avgoverfreq = 'yes';
    cfg.latency     = [0.5 1.5];
    freq            = ft_selectdata(cfg,freq);
        
    % save    
    mkdir([dir_bids,'derivatives/',subj_handle,'/eeg/'])
    save([dir_bids,'derivatives/',subj_handle,'/eeg/',subj_handle,'_task-rf_eeg-freq.mat'],'freq','brainGeo')    
    
    % clean workspace
    clear subj_handle freq avg_pow std_pow
end

%% Get Group Averages
% predefine group data matrix
group_freq = cell(n_subj,2);

% load template grid
load('Y:/projects/general/fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d4mm.mat'); 
template_grid = sourcemodel; clear sourcemodel

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject handle
    subj_handle = sprintf('sub-%02.0f',subj);
   
    % load data
    load([dir_bids,'derivatives/',subj_handle,'/eeg/',subj_handle,'_task-rf_eeg-freq.mat']);
    
    % find audio/video hit/miss trials
    [audio_bool,hit_bool]   = get_splits_eeg(freq);
    
    % split video data
    cfg                     = [];
    cfg.trials              = hit_bool==1 & audio_bool==0;
    group_freq{subj,1}      = ft_freqdescriptives(cfg, freq);

    cfg                     = [];
    cfg.trials              = hit_bool==0 & audio_bool==0;
    group_freq{subj,2}      = ft_freqdescriptives(cfg, freq);
    
    % split audio data
    cfg                     = [];
    cfg.trials              = hit_bool==1 & audio_bool==1;
    group_freq{subj,3}      = ft_freqdescriptives(cfg, freq);

    cfg                     = [];
    cfg.trials              = hit_bool==0 & audio_bool==1;
    group_freq{subj,4}      = ft_freqdescriptives(cfg, freq);
    
    % add source details for each conditon
    for condition = 1 : size(group_freq,2)
        
        % add source fields
        group_freq{subj,condition}.pos          = template_grid.pos;
        group_freq{subj,condition}.dim          = brainGeo.dim;
        group_freq{subj,condition}.inside       = brainGeo.inside;
        group_freq{subj,condition}.powdimord	= 'pos';
        group_freq{subj,condition}.pow          = group_freq{subj,condition}.powspctrm;
        group_freq{subj,condition}              = rmfield(group_freq{subj,condition},{'dimord','powspctrm','time','freq','cumtapcnt'});
        
        % add 'outside' co-ordinates to powspctrm
        x                                       = nan(size(group_freq{subj,condition}.inside));
        x(group_freq{subj,condition}.inside==1) = group_freq{subj,condition}.pow;
        group_freq{subj,condition}.pow          = x;        
    end
end
   
% get grand average
cfg                 = [];
cfg.keepindividual  = 'yes';
grand_hits.video    = ft_sourcegrandaverage(cfg, group_freq{:,1});
grand_misses.video  = ft_sourcegrandaverage(cfg, group_freq{:,2});
grand_hits.audio    = ft_sourcegrandaverage(cfg, group_freq{:,3});
grand_misses.audio  = ft_sourcegrandaverage(cfg, group_freq{:,4});

% save
mkdir([dir_bids,'derivatives/group/eeg/'])
save([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-freq.mat'],'grand_hits','grand_misses')
keep dir_root

%% Run Statistical Analysis across Whole Brain
% load data
load([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-freq.mat'])

% get modality names
modality = fieldnames(grand_hits);

% cycle through each modality
for i = 1 : numel(modality)

    % run statisitics
    cfg                   = [];
    cfg.dim               = grand_hits.(modality{i}).dim;  % specify dimensions of your source grid
    cfg.design            = [1:size(grand_hits.(modality{i}).pow,1), 1:size(grand_hits.(modality{i}).pow,1); ones(1,size(grand_hits.(modality{i}).pow,1)), ones(1,size(grand_hits.(modality{i}).pow,1))+1];
    cfg.uvar              = 1;
    cfg.ivar              = 2;
    cfg.method            = 'montecarlo';
    cfg.parameter         = 'pow';
    cfg.statistic         = 'ft_statfun_depsamplesT';
    cfg.correctm          = 'cluster';
    cfg.clusteralpha      = 0.05;
    cfg.numrandomization  = 2000;
    cfg.alpha             = 0.05;
    cfg.tail              = -1;
    cfg.clustertail       = -1;
    stat.(modality{i})    = ft_sourcestatistics(cfg, grand_hits.(modality{i}), grand_misses.(modality{i}));
    
    % extract key values
    p(i,1)  = round(stat.(modality{i}).negclusters(1).prob,3); %#ok<SAGROW>
    t(i,1)  = round(stat.(modality{i}).negclusters(1).clusterstat ./ sum(stat.(modality{i}).negclusterslabelmat == 1),3); %#ok<SAGROW>
    ci(i,1) = round(stat.(modality{i}).negclusters(1).cirange,3); %#ok<SAGROW>
    
    % calculate cohen's dz
    dz(i,1) = round(t(i,1)./ sqrt(size(grand_hits.(modality{i}).pow,1)),3); %#ok<SAGROW>
    stat.(modality{i}).cohens_dz = dz(i,1);
end

% create results table
results_table = table(modality,p,ci,t,dz) %#ok<NOPTS>

% save
save([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-stat.mat'],'stat','results_table')

%% Create Masked Source and Export
% load data
load([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-stat.mat'])

% load template MRI
mri = ft_read_mri('Y:/projects/general/fieldtrip-20170319/template/anatomy/single_subj_T1_1mm.nii');

% get modality names
modality = fieldnames(stat);

% cycle through each modality
for i = 1 : numel(modality)

    % get indices of clusters
    clus_idx = stat.(modality{i}).negclusterslabelmat==1;

    % create source data structure
    source                  = [];
    source.inside           = stat.(modality{i}).inside;
    source.dim              = stat.(modality{i}).dim;
    source.pos              = stat.(modality{i}).pos*10;
    source.unit             = 'mm';
    
    % define powspctrm of cluster
    source.pow              = nan(size(stat.(modality{i}).pos,1),1);     
    source.pow(clus_idx)	= stat.(modality{i}).stat(clus_idx); 
    
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
    cfg.filename      = [dir_bids,'data/eeg/source/stat_',modality{i},'.nii'];  % enter the desired file name
    cfg.filetype      = 'nifti';
    cfg.coordsys      = 'spm';
    cfg.vmpversion    = 2;
    cfg.datatype      = 'float';
    ft_volumewrite(cfg,source);      % be sure to use your interpolated source data
    
    % reslice to 1mm isometric to match template MRI
    reslice_nii([dir_bids,'data/eeg/source/stat_',modality{i},'.nii'],[dir_bids,'data/eeg/source/stat_',modality{i},'_sm.nii'],[1 1 1]);
end

%% Extract Raw Power of Cluster 
% load data
load([dir_bids,'data/eeg/source/grand_freq.mat'])

% get indices of clusters
clus_idx = stat.video.negclusterslabelmat==1;

% get mean magnitude of cluster
hits = nanmean(grand_hits.video.pow(:,clus_idx),2);
misses = nanmean(grand_misses.video.pow(:,clus_idx),2);
rse = hits - misses;

% create table
tbl = table(rse);

% write table
writetable(tbl,[dir_repos,'data/sourceEEG.csv'],'Delimiter',',')

%% Get Time-Series of Cluster
% load cluster data
load([dir_bids,'data/eeg/source/stat.mat'])

% predefine group data matrix
group_freq = cell(n_subj,2);

% load template grid
load('Y:/projects/general/fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d4mm.mat'); 
template_grid = sourcemodel; clear sourcemodel

% cycle through each subject
for subj = 1 : n_subj
    
    % load
    load([dir_bids,'data/eeg/source/subj',sprintf('%02.0f',subj),'_source.mat'])
    
    % find audio/video hit/miss trials
    [audio_bool,hit_bool]   = get_splits_eeg(source);

    % define electrodes in cluster
    labels = source.label(stat.video.negclusterslabelmat(stat.video.inside)==1);
    
    % get time-frequency representation of data
    cfg             = [];
    cfg.trials      = audio_bool == 0;
    cfg.channel     = labels;
    cfg.keeptrials  = 'yes';
    cfg.method      = 'wavelet';
    cfg.width       = 6;
    cfg.output      = 'pow';
    cfg.toi         = -0.5:0.05:2;
    cfg.foi         = 8:25; % 50hz sampling
    cfg.pad         = 'nextpow2';
    freq            = ft_freqanalysis(cfg, source); clear source

    % convert powspctrm to single
    freq.powspctrm  = single(freq.powspctrm);
    
    % z-transform
    avg_pow         = repmat(nanmean(nanmean(freq.powspctrm,4),1),[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
    std_pow         = repmat(nanstd(nanmean(freq.powspctrm,4),[],1),[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
    freq.powspctrm  = (freq.powspctrm - avg_pow) ./ std_pow; clear avg_pow std_pow
    
    % apply gaussian smoothing
    cfg             = [];
    cfg.fwhm_t      = 0.2;
    cfg.fwhm_f      = 2;
    freq            = smooth_TF_GA(cfg,freq);
    
    % average over time and frequency
    cfg             = [];
    cfg.avgoverfreq = 'yes';
    cfg.avgoverchan = 'yes';
    freq            = ft_selectdata(cfg,freq);
    
    % split video data
    hits_series(subj,:) = squeeze(mean(freq.powspctrm(hit_bool(audio_bool==0)==1,:,:,:),1));
    miss_series(subj,:) = squeeze(mean(freq.powspctrm(hit_bool(audio_bool==0)==0,:,:,:),1));
end

% save
csvwrite([dir_repos,'data/sourceEEG_hitSeries.csv'],hits_series)
csvwrite([dir_repos,'data/sourceEEG_missSeries.csv'],miss_series)
   
