%% EEG: Source Analysis
% prepare workspace
clearvars
close all
clc

% define root directory
if ispc;    dir_root = 'Y:/projects/reinstatement_fidelity/';
            dir_tool = 'Y:/projects/general';
            dir_repos  = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
else;       dir_root = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity';
            dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general';
end

% add subfunctions
addpath([dir_root,'scripts/subfunctions'])

% define number of subjects
n_subj = 21;

%% Create AAL Source Grids
% load template grid
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d4mm.mat']); 
template_grid = sourcemodel; clear sourcemodel % rename grid

% load AAL atlas
mri = ft_read_mri([dir_tool,'general/fieldtrip-20170319/template/atlas/aal/ROI_MNI_V4.nii']);
 
% change MRI to binary 'in-atlas' vs. 'out-atlas'
mri.anatomy = double(mri.anatomy > 0); 

% interpolate mri with grid
cfg                 = [];
cfg.parameter       = 'anatomy';
cfg.interpmethod    = 'nearest';
roi                 = ft_sourceinterpolate(cfg,mri,template_grid);

% clean up
clear mri template_grid cfg

%% Decompose Source Signal
% cycle through each subject
for subj = 1 : n_subj
        
    % load data
    load([dir_root,'data/eeg/preprocessing/subj',sprintf('%02.0f',subj),'_data.mat'])
        
    % load electrode layout, grid and boundary element model    
    load([dir_root,'data/eeg/headmodels/subj',sprintf('%02.0f',subj),'/elec.mat'])
    load([dir_root,'data/eeg/headmodels/subj',sprintf('%02.0f',subj),'/vol.mat'])
    load([dir_root,'data/eeg/headmodels/subj',sprintf('%02.0f',subj),'/leadfield.mat'])
    
    % restrict grid to AAL atlas
    grid.inside = roi.anatomy(:) > 0;
    
    % downsample data
    cfg             = [];
    cfg.resamplefs  = 50; % save some space (limits later analyses to 25Hz)
    data            = ft_resampledata(cfg, data);

    % timelock data
    cfg                     = [];
    cfg.covariancewindow    = 'all';
    cfg.covariance          = 'yes';
    cfg.keeptrials          = 'yes';
    tml                     = ft_timelockanalysis(cfg, data);

    % get filters for video and audio conditions
    cfg                 = [];
    cfg.method          = 'lcmv';
    cfg.lcmv.keepfilter = 'yes';
    cfg.lcmv.fixedori   = 'yes';
    cfg.lcmv.lambda     = '5%';
    cfg.headmodel       = vol;
    cfg.elec            = elec;
    cfg.grid            = grid;
    tmp                 = ft_sourceanalysis(cfg, tml);
    filter              = cell2mat(tmp.avg.filter(tmp.inside));
    
    % store brain geometry
    brainGeo.dim        = tmp.dim;
    brainGeo.inside     = tmp.inside;
    brainGeo.pos        = tmp.pos;

    % combine filters with raw video data to "beam" scalp video data to source space
    source            = [];
    source.time       = data.time;       % transfer time information
    source.trialinfo  = data.trialinfo;  % transfer trial information    
    for c = 1 : sum(grid.inside); source.label{c,1} = ['S' num2str(c)]; end % create labels for each virtual electrode    
    for j = 1 : numel(data.trial); source.trial{1,j} = single(filter*data.trial{1,j}); end % for each trial, apply filters to the recorded data
    
    % save
    save([dir_root,'data/eeg/source/subj',sprintf('%02.0f',subj),'_source.mat'],'source','brainGeo','-v7.3')
    clear cfg aal data tml tmp filter template_grid elec grid vol brainGeo source c j
    
end

%% Get TF of Source Data
% predefine group data matrix
group_freq = cell(n_subj,2);

% load template grid
load('Y:/projects/general/fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d4mm.mat'); 
template_grid = sourcemodel; clear sourcemodel

% cycle through each subject
for subj = 1 : n_subj
    
    % load
    load([dir_root,'data/eeg/source/subj',sprintf('%02.0f',subj),'_source.mat'])
    
    % find audio/video hit/miss trials
    [audio_bool,hit_bool]   = get_splits_eeg(source);

    % get time-frequency representation of data
    cfg             = [];
    cfg.keeptrials  = 'yes';
    cfg.method      = 'wavelet';
    cfg.width       = 6;
    cfg.output      = 'pow';
    cfg.toi         = -1:0.25:3;
    cfg.foi         = 8:2:24; % 50hz sampling
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
    save([dir_root,'data/eeg/source/subj',sprintf('%02.0f',subj),'_freq.mat'],'freq')
    
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
    
    % save
    save([dir_root,'data/eeg/source/subj',sprintf('%02.0f',subj),'_freq_conditions.mat'],'freq')
    clear cfg freq freqtmp x brainGeo audio_bool hit_bool condition
end
   
% get grand average
cfg                 = [];
cfg.keepindividual  = 'yes';
grand_hits.video    = ft_sourcegrandaverage(cfg, group_freq{:,1});
grand_misses.video  = ft_sourcegrandaverage(cfg, group_freq{:,2});
grand_hits.audio    = ft_sourcegrandaverage(cfg, group_freq{:,3});
grand_misses.audio  = ft_sourcegrandaverage(cfg, group_freq{:,4});

% save
save([dir_root,'data/eeg/source/grand_freq.mat'],'grand_hits','grand_misses')
keep dir_root

%% Run Statistical Analysis across Whole Brain
% load data
load([dir_root,'data/eeg/source/grand_freq.mat'])

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
save([dir_root,'data/eeg/source/stat.mat'],'stat')

%% Create Masked Source and Export
% load data
load([dir_root,'data/eeg/source/stat.mat'])

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
    cfg.filename      = [dir_root,'data/eeg/source/stat_',modality{i},'.nii'];  % enter the desired file name
    cfg.filetype      = 'nifti';
    cfg.coordsys      = 'spm';
    cfg.vmpversion    = 2;
    cfg.datatype      = 'float';
    ft_volumewrite(cfg,source);      % be sure to use your interpolated source data
    
    % reslice to 1mm isometric to match template MRI
    reslice_nii([dir_root,'data/eeg/source/stat_',modality{i},'.nii'],[dir_root,'data/eeg/source/stat_',modality{i},'_sm.nii'],[1 1 1]);
end

%% Extract Raw Power of Cluster 
% load data
load([dir_root,'data/eeg/source/grand_freq.mat'])

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
load([dir_root,'data/eeg/source/stat.mat'])

% predefine group data matrix
group_freq = cell(n_subj,2);

% load template grid
load('Y:/projects/general/fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d4mm.mat'); 
template_grid = sourcemodel; clear sourcemodel

% cycle through each subject
for subj = 1 : n_subj
    
    % load
    load([dir_root,'data/eeg/source/subj',sprintf('%02.0f',subj),'_source.mat'])
    
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
   
