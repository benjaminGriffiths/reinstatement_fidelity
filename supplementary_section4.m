%% EEG: Source Analysis
% prepare workspace
clearvars
close all
clc

% define root directory
if ispc;    dir_bids = 'Y:/projects/reinstatement_fidelity/bids_data/';
            dir_tool = 'Y:/projects/general/';
            dir_repos = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
else        dir_bids = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/bids_data/';
            dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';
            dir_repos = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/scripts/';
end

% add subfunctions
addpath([dir_repos,'subfunctions'])

% define number of subjects
n_subj = 21;

%% Get TF of Source Data
% load stat
load([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-stat.mat'])

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject handle
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_bids,'sourcedata/',subj_handle,'/eeg/'];
    
    % load
    load([dir_bids,'sourcedata/',subj_handle,'/eeg/',subj_handle,'_task-rf_eeg.mat'])
        
    % find audio/video hit/miss trials
    [audio_bool,hit_bool]   = get_splits_eeg(source);
        
    % restrict source to channels in cluster
    cfg             = [];
    cfg.trials      = find(audio_bool == 0 & hit_bool == 1);
    cfg.channel     = source.label(stat.video.negclusterslabelmat(stat.video.inside==1)==1);
    source          = ft_selectdata(cfg,source);
    
    % get time-frequency representation of data
    cfg             = [];
    cfg.keeptrials  = 'yes';
    cfg.method      = 'wavelet';
    cfg.width       = 6;
    cfg.output      = 'pow';
    cfg.toi         = 0.5:0.25:1.5;
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
    freq            = ft_selectdata(cfg,freq);
     
    % extract confidence values
    for trl = 1 : numel(freq.trialinfo)
        confi(trl,1) = freq.trialinfo{trl}.confidenceScore;
    end
        
    % cycle through each channel and correlate
    for i = 1 : size(freq.powspctrm,2)
        
        % extract ranked power
        X = freq.powspctrm(:,i);
        
        % chuck out nans
        X = X(~isnan(confi));
        Y = confi(~isnan(confi));
        
        % get spearmans correlation
        r(i,1) = corr(X,Y,'type','Spearman');
    end
    
    % add to freq structure
    freq = rmfield(freq,{'cumtapcnt','trialinfo'});
    freq.powspctrm = r;
    freq.dimord = 'chan_freq_time';
    
    % save    
    mkdir([dir_bids,'derivatives/',subj_handle,'/eeg/'])
    save([dir_bids,'derivatives/',subj_handle,'/eeg/',subj_handle,'_task-rf_eeg-confidence.mat'],'freq','brainGeo')    
    
    % clean workspace
    clear subj_handle freq cfg i X Y r confi trl
end

%% Get Group Averages
% load stat
load([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-stat.mat'])

% predefine group data matrix
group_freq = cell(n_subj,1);

% load template grid
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d4mm.mat']); 
template_grid = sourcemodel; clear sourcemodel

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject handle
    subj_handle = sprintf('sub-%02.0f',subj);
   
    % load data
    load([dir_bids,'derivatives/',subj_handle,'/eeg/',subj_handle,'_task-rf_eeg-confidence.mat']);

    % add data to group
    group_freq{subj,1} = freq;
    
    % add source fields
    group_freq{subj,1}.pos          = template_grid.pos;
    group_freq{subj,1}.dim          = brainGeo.dim;
    group_freq{subj,1}.inside       = stat.video.negclusterslabelmat==1;
    group_freq{subj,1}.powdimord	= 'pos';
    group_freq{subj,1}.pow          = group_freq{subj,1}.powspctrm;
    group_freq{subj,1}              = rmfield(group_freq{subj,1},{'dimord','powspctrm','time','freq'});

    % add 'outside' co-ordinates to powspctrm
    x                               = nan(size(group_freq{subj,1}.inside));
    x(group_freq{subj,1}.inside==1) = group_freq{subj,1}.pow;
    group_freq{subj,1}.pow          = x;  
    
    % tidy up
    clear x freq subj_handle
end
   
% get grand average
cfg                 = [];
cfg.keepindividual  = 'yes';
grand_freq          = ft_sourcegrandaverage(cfg, group_freq{:,1});

% save
mkdir([dir_bids,'derivatives/group/eeg/'])
save([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-confidence.mat'],'grand_freq')
keep dir_root dir_bids dir_tool

%% Run Statistical Analysis across Whole Brain
% load data
load([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-confidence.mat'])

% create null hypothesis
grand_null = grand_freq;
grand_null.pow = zeros(size(grand_freq.pow));

% run statisitics
cfg                   = [];
cfg.dim               = grand_freq.dim;  % specify dimensions of your source grid
cfg.design            = [1:size(grand_freq.pow,1), 1:size(grand_freq.pow,1); ones(1,size(grand_freq.pow,1)), ones(1,size(grand_freq.pow,1))+1];
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
stat                  = ft_sourcestatistics(cfg, grand_freq, grand_null);

% extract key values
p  = round(stat.negclusters(1).prob,3);
t  = round(stat.negclusters(1).clusterstat ./ sum(stat.negclusterslabelmat == 1),3);
ci = round(stat.negclusters(1).cirange,3);

% calculate cohen's dz
dz = round(t./ sqrt(size(grand_freq.pow,1)),3);
stat.cohens_dz = dz;

% create results table
results_table = table(p,ci,t,dz) %#ok<NOPTS>

% save
save([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-confistat.mat'],'stat','results_table')


