%% fMRI Univariate Analysis
% prepare workspace
clearvars
close all
clc

% initalise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% define directories
git_dir  = 'E:\bjg335\projects\reinstatement_fidelity\';
bear_dir = 'Y:\projects\reinstatement_fidelity\bids_data\';
data_dir = 'E:\bids\';

% copy bids data to local directory (faster computation)
copyfile(bear_dir,data_dir)

% add subfunctions
addpath([git_dir,'subfunctions'])

% define key parameters
n_runs      = 8;
n_volumes   = 255;
n_slices    = 32;
TR          = 2;
n_subj      = 21;

% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    subj_dir = [data_dir,subj_handle,'\'];
    