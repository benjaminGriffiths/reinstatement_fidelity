%% Power Analysis
% clear workspace
clearvars
close all
clc

% define root directory
if ispc;    dir_bids = 'Y:/projects/reinstatement_fidelity/bids_data/';
            dir_tool = 'Y:/projects/general/';
            dir_repos = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
else;       dir_bids = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/bids_data/';
            dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';
            dir_repos = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/git_clone/reinstatement_fidelity/'; % repository directory
end

% add subfunctions
addpath([dir_repos,'subfunctions'])

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

%% Get Timelocked Representation of Each Condition
% predefine cell for group data
group_freq = cell(n_subj,2);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];
    
    % load in raw data
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'])  
    
    % get time-frequency representaiton of data for specified conditions
    group_freq{subj,1} = get_basediff_timefrequency(source,'encoding','visual');
    group_freq{subj,2} = get_memdiff_timefrequency(source,'retrieval','visual');
end
    
% predefine cells to house concatenated group data
grand_freq = cell(size(group_freq,2),1);

% cycle through conditions
for i = 1 : size(group_freq,2)
    
    % get grand average of subjects
    cfg                 = [];
    cfg.keepindividual  = 'yes';
    grand_freq{i,1}     = ft_freqgrandaverage(cfg,group_freq{:,i});
    
    % collapse over time/frequency
    grand_freq{i,1}.pow = mean(mean(grand_freq{i,1}.powspctrm,4),3);
    
    % add in source model details
    grand_freq{i,1}.powdimord   = 'rpt_pos';
    grand_freq{i,1}.dim         = roi.dim;
    grand_freq{i,1}.inside      = roi.inside(:);
    grand_freq{i,1}.pos         = roi.pos;
    grand_freq{i,1}.cfg         = [];
    
    % add in "outside" virtual electrods
    tmp_pow = zeros(size(grand_freq{i,1}.pow,1),size(grand_freq{i,1}.inside,1));
    tmp_pow(:,grand_freq{i,1}.inside) = grand_freq{i,1}.pow;
    grand_freq{i,1}.pow = tmp_pow;
    
    % remove freq details
    grand_freq{i,1} = rmfield(grand_freq{i,1},{'powspctrm','label','freq','time','dimord'});
end
    
% save data
mkdir([dir_bids,'derivatives/group/eeg/'])
save([dir_bids,'derivatives/group/eeg/group_task-all_eeg-source.mat'],'grand_freq');

%% Run Statistics
% predefine cell for statistics
cfg.tail        = -1;
cfg.parameter   = 'pow';
[stat,tbl]      = run_oneSampleT(cfg, grand_freq);

% save data
save([dir_bids,'derivatives/group/eeg/group_task-all_eeg-stat.mat'],'stat','tbl');

%% Extract Raw Power of Cluster 
% define cluster names
plot_names = {'encoding','retrieval'};

% prepare table for stat values
tbl = array2table(zeros(n_subj,2),'VariableNames',plot_names);

% cycle through conditions
for i = 1 : numel(stat)
     
    % get indices of clusters
    clus_idx = stat{i}.negclusterslabelmat==1;

    % create table
    tbl.(plot_names{i})(:,1) = nanmean(grand_freq{i}.pow(:,clus_idx),2);
end

% write table
writetable(tbl,[dir_bids,'derivatives/group/eeg/group_task-all_eeg-cluster.csv'],'Delimiter',',')
copyfile([dir_bids,'derivatives/group/eeg/group_task-all_eeg-cluster.csv'],[dir_repos,'data/fig2_data/group_task-rf_eeg-cluster.csv'])

%% Create Surface Plots
% load template MRI
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat']);

% define plot names
plot_names = {'encoding','retrieval'};

% cycle through conditions
for i = 1 : numel(stat)
        
    % get indices of clusters
    clus_idx = stat{i}.negclusterslabelmat==1;

    % create source data structure
    source                  = [];
    source.inside           = stat{i}.inside;
    source.dim              = stat{i}.dim;
    source.pos              = stat{i}.pos*10;
    source.unit             = 'mm';
    
    % define powspctrm of cluster
    source.pow              = nan(size(stat{i}.pos,1),1);     
    source.pow(clus_idx)	= stat{i}.stat(clus_idx); 
    
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
    cfg.filename      = [dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-map.nii'];  % enter the desired file name
    cfg.filetype      = 'nifti';
    cfg.coordsys      = 'spm';
    cfg.vmpversion    = 2;
    cfg.datatype      = 'float';
    ft_volumewrite(cfg,source);      % be sure to use your interpolated source data
    
    % reslice to 1mm isometric to match template MRI
    reslice_nii([dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-map.nii'],[dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-map.nii'],[1 1 1]);
    copyfile([dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-map.nii'],[dir_repos,'data/fig2_data/group_task-',plot_names{i},'_eeg-map.nii'])    
end

%% Collect TF of Cluster
% cycle through each subject
for subj = 1 : n_subj 
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];
    
    % load in raw data
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'])  
    data = source; clear source
    
    % predefine arrays for memory performance
    mem_performance = zeros(numel(data.trial),1);
    operation_to_include = zeros(numel(data.trial),1);
    modality_to_include = zeros(numel(data.trial),1);
    
    % cycle through each trial
    for trl = 1 : numel(data.trial)

        % mark trials that do match specified operation
        operation_to_include(trl) = strcmpi(data.trialinfo{trl}.operation,'encoding');
        
        % mark trials that do match specified operation
        modality_to_include(trl) = strcmpi(data.trialinfo{trl}.modality,'visual');

        % get memory
        mem_performance(trl) = data.trialinfo{trl}.recalled;
        
    end

    % get time-frequency
    cfg                 = [];
    cfg.keeptrials      = 'yes';
    cfg.method          = 'wavelet';
    cfg.width           = 6;
    cfg.output          = 'pow';
    cfg.toi             = -2:0.05:2;
    cfg.foi             = 3:0.5:40;
    cfg.pad             = 'nextpow2';
    
    % get encoding data
    cfg.channel         = data.label(stat{1}.negclusterslabelmat(stat{1}.inside)==1);
    cfg.trials          = find(operation_to_include == 1 & modality_to_include == 1);
    freq                = ft_freqanalysis(cfg, data);
    
    % get time-averaged mean and standard deviation of power for each channel and frequency
    raw_pow = mean(freq.powspctrm,4);
    avg_pow = mean(raw_pow,1);
    std_pow = std(raw_pow,[],1);

    % replicate matrices to match freq.powspctrm
    avg_pow = repmat(avg_pow,[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
    std_pow = repmat(std_pow,[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
    
    % z-transform power
    freq.powspctrm = (freq.powspctrm - avg_pow) ./ std_pow;
    clear raw_pow avg_pow std_pow
    
    % smooth
    config = [];
    config.fwhm_t=0.2;
    config.fwhm_f=2;
    freq = smooth_TF_GA(config,freq);

    % get pre-post difference
    percept_diff(subj,:,:) = squeeze(mean(mean(freq.powspctrm,2),1));
    
    % get retrieval data
    cfg.channel         = data.label(stat{2}.negclusterslabelmat(stat{1}.inside)==1);
    cfg.trials          = find(operation_to_include == 0 & modality_to_include == 1);
    freq                = ft_freqanalysis(cfg, data);
    
    % get time-averaged mean and standard deviation of power for each channel and frequency
    raw_pow = mean(freq.powspctrm,4);
    avg_pow = mean(raw_pow,1);
    std_pow = std(raw_pow,[],1);

    % replicate matrices to match freq.powspctrm
    avg_pow = repmat(avg_pow,[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
    std_pow = repmat(std_pow,[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
    
    % z-transform power
    freq.powspctrm = (freq.powspctrm - avg_pow) ./ std_pow;
    clear raw_pow avg_pow std_pow

    % smooth
    cfg = [];
    cfg.fwhm_t=0.2;
    cfg.fwhm_f=2;
    freq = smooth_TF_GA(cfg,freq);
    
    % get memory difference
    vismem = mem_performance(operation_to_include == 0 & modality_to_include == 1);
    memory_diff(subj,:,:,1) = squeeze(mean(mean(freq.powspctrm(vismem==1,:,:,:),1),2)); 
    memory_diff(subj,:,:,2) = squeeze(mean(mean(freq.powspctrm(vismem==0,:,:,:),1),2));
end

% extract variables of interest
percept_freq(:,:,1) = squeeze(mean(percept_diff(:,:,freq.time>=0.5 & freq.time<=1.5),3));
percept_freq(:,:,2) = squeeze(mean(percept_diff(:,:,freq.time>=-1 & freq.time<=0),3));
percept_time        = squeeze(mean(percept_diff(:,freq.freq>=8 & freq.freq<=30,:),2));

% extract variables of interest
memory_freq  = squeeze(mean(memory_diff(:,:,freq.time>=0.5 & freq.time<=1.5,:),3));
memory_time  = squeeze(mean(memory_diff(:,freq.freq>=8 & freq.freq<=30,:,:),2));

% reshape memory variables
memory_freq = reshape(memory_freq,[21 150]);
memory_time = reshape(memory_time,[21 122]);
percept_freq = reshape(percept_freq,[21 150]);

% save as csv
csvwrite([dir_repos,'data/fig2_data/group_task-percept_eeg-freqseries.csv'],percept_freq)
csvwrite([dir_repos,'data/fig2_data/group_task-percept_eeg-timeseries.csv'],percept_time)
csvwrite([dir_repos,'data/fig2_data/group_task-memory_eeg-freqseries.csv'],memory_freq)
csvwrite([dir_repos,'data/fig2_data/group_task-memory_eeg-timeseries.csv'],memory_time)
