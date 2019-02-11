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

%% Prepare Occipital ROI
% load template grid
load('Y:/projects/general/fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat'); 
template_grid = sourcemodel;

% load whole brain atlas
mri = ft_read_mri([dir_bids,'sourcedata/masks/whole_brain.nii']);
 
% interpolate clusters with grid
cfg             = [];
cfg.parameter	= 'anatomy';
roi             = ft_sourceinterpolate(cfg,mri,sourcemodel);

% define additional roi parameters
roi.inside      = ~isnan(roi.anatomy) & roi.anatomy > 0;
roi.anatomy     = double(~isnan(roi.anatomy) & roi.anatomy > 0);
roi.pos         = sourcemodel.pos;

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
roi2                = ft_sourceinterpolate(cfg,mri,template_grid);

% determine inside-inside
roi.insideRoi = roi2.anatomy(roi.inside(:)==1)>0;

% clean up
clear mri cfg

%% Get Timelocked Representation of Each Condition
% predefine cell for group data
group_erp = zeros(n_subj,121);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];
    
    % load in raw data
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'])  
    data = source; clear source
    
    % predefine conditional arrays to include all trials
    operation_to_include = zeros(numel(data.trial),1);
    modality_to_include  = zeros(numel(data.trial),1);

    % cycle through each trial
    for trl = 1 : numel(data.trial)

        % mark trials that do match specified operation
        operation_to_include(trl) = strcmpi(data.trialinfo{trl}.operation,'encoding');

        % mark trials that do match specified operation
        modality_to_include(trl) = strcmpi(data.trialinfo{trl}.modality,'visual');
    end

    % select data
    cfg                 = [];
    cfg.channel         = data.label(roi.insideRoi==1);
    cfg.trials          = operation_to_include == 1 & modality_to_include == 1;
    cfg.latency         = [-0.2 1];
    occip_data          = ft_selectdata(cfg,data);

    % preprocess
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [8 30];
    data = ft_preprocessing(cfg,data);
    
    % timelock data
    cfg = [];
    tml = ft_timelockanalysis(cfg,occip_data);
    
    % baseline correct
    cfg = [];
    cfg.baseline = [-0.2 0];
    tml = ft_timelockbaseline(cfg,tml);
    
    % get group erp
    group_erp(subj,:) = mean(tml.avg,1);
end
