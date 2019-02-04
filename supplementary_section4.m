%% EEG: Source Analysis
% prepare workspace
% initialisation
clearvars
clc

% define root directory
if ispc        
    dir_root = 'Y:/projects/reinstatement_fidelity/';
    dir_bids = [dir_root,'bids_data/'];
    dir_tool = 'Y:/projects/general/';
    dir_repos = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
    
    % add subfunctions
    addpath([dir_repos,'subfunctions'])
    
else       
    dir_root = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/';
    dir_bids = [dir_root,'bids_data/'];
    dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';            
    dir_repos = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/git_clone/reinstatement_fidelity/'; % repository directory
    
    % add subfunctions
    addpath([dir_repos,'subfunctions'])
end

% define number of subjects
n_subj = 21;

%% Prepare ROI
% load sourcemodel
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat']);

% load AAL atlas
mri = ft_read_mri([dir_bids,'sourcedata/masks/whole_brain.nii']);
 
% interpolate clusters with grid
cfg             = [];
cfg.parameter	= 'anatomy';
roi             = ft_sourceinterpolate(cfg,mri,sourcemodel);

% define additional roi parameters
roi.inside      = ~isnan(roi.anatomy) & roi.anatomy > 0;
roi.anatomy     = double(~isnan(roi.anatomy) & roi.anatomy > 0);
roi.pos         = sourcemodel.pos;

%% Get TF of Source Data
% cycle through each subject
for subj = 1 : n_subj
    
    % define subject handle
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_bids,'sourcedata/',subj_handle,'/eeg/'];
    
    % load
    load([dir_bids,'sourcedata/',subj_handle,'/eeg/',subj_handle,'_task-rf_eeg-source.mat'])
        
    % cycle through each trial
    isVisRet = nan(size(source.trial));
    for trl = 1 : numel(source.trial)
        isVisRet(trl) = strcmpi(source.trialinfo{trl}.operation,'retrieval')&strcmpi(source.trialinfo{trl}.modality,'visual');
    end
       
    % restrict source to channels in cluster
    cfg             = [];
    cfg.trials      = find(isVisRet);
    source          = ft_selectdata(cfg,source);
    
    % get time-frequency representation of data
    cfg             = [];
    cfg.keeptrials  = 'yes';
    cfg.method      = 'wavelet';
    cfg.width       = 6;
    cfg.output      = 'pow';
    cfg.toi         = 0.5:0.1:1.5;
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
    confi = nan(size(freq.trialinfo));
    for trl = 1 : numel(freq.trialinfo)
        confi(trl,1) = freq.trialinfo{trl}.confidence;
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
    save([dir_bids,'derivatives/',subj_handle,'/eeg/',subj_handle,'_task-rf_eeg-confidence.mat'],'freq')    
    
    % clean workspace
    clear subj_handle freq cfg i X Y r confi trl
end

%% Get Group Averages
% predefine group data matrix
group_freq = cell(n_subj,1);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject handle
    subj_handle = sprintf('sub-%02.0f',subj);
   
    % load data
    load([dir_bids,'derivatives/',subj_handle,'/eeg/',subj_handle,'_task-rf_eeg-confidence.mat']);

    % add data to group
    group_freq{subj,1} = freq;
    
    % add source fields
    group_freq{subj,1}.pos          = roi.pos;
    group_freq{subj,1}.dim          = roi.dim;
    group_freq{subj,1}.inside       = roi.inside(:);
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
save([dir_bids,'derivatives/group/eeg/group_task-confidence_eeg-confidence.mat'],'grand_freq')

%% Run Statistics
% predefine cell for statistics
cfg.tail        = -1;
cfg.parameter   = 'pow';
[stat,tbl]      = run_oneSampleT(cfg, grand_freq);

% save data
save([dir_bids,'derivatives/group/eeg/group_task-confidence_eeg-stat.mat'],'stat','tbl');


