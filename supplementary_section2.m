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

% load whole brain
mri = ft_read_mri([dir_bids,'sourcedata/masks/whole_brain.nii']);

% interpolate clusters with grid
cfg             = [];
cfg.parameter	= 'anatomy';
roi             = ft_sourceinterpolate(cfg,mri,sourcemodel);

% define additional roi parameters
roi.inside      = roi.anatomy(:) > 0;
roi.pos         = sourcemodel.pos;

%% Get Timelocked Representation of Each Condition
% predefine cell for group data
group_erp = zeros(n_subj,1817,121);

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
    cfg.trials          = operation_to_include == 1 & modality_to_include == 1;
    cfg.latency         = [-0.2 1];
    data                = ft_selectdata(cfg,data);

    % timelock data
    cfg                 = [];
    cfg.keeptrials      = 'yes';
    tml                 = ft_timelockanalysis(cfg,data);
    
    % baseline correct
    tml.trial = tml.trial - repmat(mean(tml.trial(:,:,tml.time<0),3),[1 1 size(tml.trial,3)]);
    
    % get group erp
    group_erp(subj,:,:) = mean(tml.trial,1);
end

%% Define Occipital ROI
% load AAL
mri = ft_read_mri([dir_bids,'sourcedata/masks/AAL.nii']);

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
occipital_roi       = ft_sourceinterpolate(cfg,mri,sourcemodel);

% determine inside
occipital_roi.anatomy = occipital_roi.anatomy(:);
occipital_roi.inside  = occipital_roi.anatomy>0 & roi.inside>0;

% get labels inside
coi = occipital_roi.anatomy(roi.inside>0)==1;

%% Get ERP
time = tml.time;
avg  = squeeze(mean(mean(group_erp(:,coi,:),2),1));
plot(time,avg);


