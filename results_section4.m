%% Power Similarity Analysis
% clear workspace
clearvars
close all
clc

% define root directory
if ispc;    dir_root = 'Y:/projects/reinstatement_fidelity/';
            dir_bids = 'Y:/projects/reinstatement_fidelity/bids_data/';
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

% define mask names
mask_names = {'percept','ers'};

%% Prepare ROI
% load sourcemodel
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat']);

% load AAL atlas
mri_aal = ft_read_mri([dir_root,'bids_data/sourcedata/masks/whole_brain.nii']);

% interpolate clusters with grid
cfg             = [];
cfg.parameter	= 'anatomy';
wb_roi          = ft_sourceinterpolate(cfg,mri_aal,sourcemodel);

% cycle through masks
for mask = 1 : numel(mask_names)

    % load mask
    mri_mask = ft_read_mri([dir_root,'bids_data/derivatives/group/rsa-',mask_names{mask},'/grand_cluster_dilated.nii']);

    % interpolate clusters with grid
    cfg             = [];
    cfg.parameter	= 'anatomy';
    roi{mask}       = ft_sourceinterpolate(cfg,mri_mask,sourcemodel);

    % define additional roi parameters
    roi{mask}.inside      = ~isnan(roi{mask}.anatomy) & roi{mask}.anatomy > 0;
    roi{mask}.insideWB    = ~isnan(wb_roi.anatomy) & wb_roi.anatomy > 0;
    roi{mask}.anatomy     = double(~isnan(roi{mask}.anatomy) & roi{mask}.anatomy > 0);
    roi{mask}.pos         = sourcemodel.pos;
end

%% Get Timelocked Representation of Each Condition
% predefine cell for group data
group_freq = cell(n_subj,1);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];
    
    % load in raw data
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'])  
    
    % restrict to roi
    enc_source = remove_source_electrodes(source,roi{1});
    ret_source = remove_source_electrodes(source,roi{2});
    
    % get time-frequency representaiton of data for specified conditions
    group_freq{subj,1} = get_powersim(enc_source,'encoding','visual');
    group_freq{subj,2} = get_powersim(ret_source,'retrieval','visual');
    
    % update command line
    fprintf('Subject %02.0f of %02.0f complete...\n',subj,n_subj)
end
    
% predefine cells to house concatenated group data
grand_freq = cell(size(group_freq,2),1);

% cycle through conditions
for i = 1 : size(group_freq,2)
    
    % get grand average of subjects
    cfg                 = [];
    cfg.keepindividual  = 'yes';
    grand_freq{i,1}     = ft_freqgrandaverage(cfg,group_freq{:,i});
    
end
    
% save data
mkdir([dir_bids,'derivatives/group/eeg/'])
save([dir_bids,'derivatives/group/eeg/group_task-all_eeg-powersim.mat'],'grand_freq');

%% Run Statistics
% predefine cell for statistics
cfg.tail            = 1;
cfg.parameter       = 'powspctrm';
[stat,tbl]          = run_oneSampleT(cfg, grand_freq);

% save data
save([dir_bids,'derivatives/group/eeg/group_task-all_eeg-powerstat.mat'],'stat','tbl');
