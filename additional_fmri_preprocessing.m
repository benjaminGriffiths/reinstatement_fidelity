%% fMRI Preprocessing
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

% get all 'sub-xx' folders in BIDS directory
files_to_copy = dir([bear_dir,'sub*']);

% copy bids data to local directory (faster computation)
for i = 1 : numel(files_to_copy)
    copyfile([files_to_copy(i).folder,'\',files_to_copy(i).name],[data_dir,files_to_copy(i).name])
end

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
    
    % slice timing correction
    matlabbatch{1}.spm.temporal.st.nslices                          = n_slices;
    matlabbatch{1}.spm.temporal.st.tr                               = TR;
    matlabbatch{1}.spm.temporal.st.ta                               = TR - (n_slices/2);
    matlabbatch{1}.spm.temporal.st.so                               = n_slices : -1 : 1;
    matlabbatch{1}.spm.temporal.st.refslice                         = n_slices / 2;
    matlabbatch{1}.spm.temporal.st.prefix                           = 'a';
    matlabbatch{1}.spm.temporal.st.scans                            = get_functional_files(subj_dir,'sub');
    
    % realign and unwarp
    matlabbatch{2}.spm.spatial.realignunwarp.data.pmscan            = '';
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.quality       = 0.9;
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.sep           = 4;
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.fwhm          = 5;
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.rtm           = 0;
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.einterp       = 2;
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.ewrap         = [0 0 0];
    matlabbatch{2}.spm.spatial.realignunwarp.eoptions.weight        = '';
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.basfcn      = [12 12];
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.regorder    = 1;
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.lambda      = 1;
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.jm          = 0;
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.fot         = [4 5];
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.sot         = [];
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.uwfwhm      = 4;
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.rem         = 1;
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.noi         = 5;
    matlabbatch{2}.spm.spatial.realignunwarp.uweoptions.expround    = 'Average';
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.uwwhich     = [2 1];
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.rinterp     = 4;
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.wrap        = [0 0 0];
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.mask        = 1;
    matlabbatch{2}.spm.spatial.realignunwarp.uwroptions.prefix      = 'u';
    matlabbatch{2}.spm.spatial.realignunwarp.data.scans             = update_functional_files(matlabbatch{1}.spm.temporal.st.scans,'a');
    
    % coregister
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.cost_fun     = 'nmi';
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.sep          = [4 2];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.tol          = [zeros(1,3)+0.02, zeros(1,9)+0.001];
    matlabbatch{3}.spm.spatial.coreg.estimate.eoptions.fwhm         = [7 7];
    matlabbatch{3}.spm.spatial.coreg.estimate.other                 = {''};  
    matlabbatch{3}.spm.spatial.coreg.estimate.ref                   = {[subj_dir,'func\meanua',subj_handle,'_task-rf_run-1_bold.nii,1']};
    matlabbatch{3}.spm.spatial.coreg.estimate.source                = {[subj_dir,'anat\',subj_handle,'_T1w.nii']};
    
    % segment
    matlabbatch{4}.spm.spatial.preproc.channel.biasreg              = 0.001;
    matlabbatch{4}.spm.spatial.preproc.channel.biasfwhm             = 60;
    matlabbatch{4}.spm.spatial.preproc.channel.write                = [0 1];
    matlabbatch{4}.spm.spatial.preproc.tissue(1).tpm                = {'C:\localToolboxes\SPM12\tpm\TPM.nii,1'};
    matlabbatch{4}.spm.spatial.preproc.tissue(2).tpm                = {'C:\localToolboxes\SPM12\tpm\TPM.nii,2'};
    matlabbatch{4}.spm.spatial.preproc.tissue(3).tpm                = {'C:\localToolboxes\SPM12\tpm\TPM.nii,3'};
    matlabbatch{4}.spm.spatial.preproc.tissue(4).tpm                = {'C:\localToolboxes\SPM12\tpm\TPM.nii,4'};
    matlabbatch{4}.spm.spatial.preproc.tissue(5).tpm                = {'C:\localToolboxes\SPM12\tpm\TPM.nii,5'};
    matlabbatch{4}.spm.spatial.preproc.tissue(6).tpm                = {'C:\localToolboxes\SPM12\tpm\TPM.nii,6'};
    matlabbatch{4}.spm.spatial.preproc.tissue(1).ngaus              = 1;
    matlabbatch{4}.spm.spatial.preproc.tissue(2).ngaus              = 1;
    matlabbatch{4}.spm.spatial.preproc.tissue(3).ngaus              = 2;
    matlabbatch{4}.spm.spatial.preproc.tissue(4).ngaus              = 3;
    matlabbatch{4}.spm.spatial.preproc.tissue(5).ngaus              = 4;
    matlabbatch{4}.spm.spatial.preproc.tissue(6).ngaus              = 2;
    matlabbatch{4}.spm.spatial.preproc.tissue(1).native             = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(2).native             = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(3).native             = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(4).native             = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(5).native             = [1 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(6).native             = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(1).warped             = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(2).warped             = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(3).warped             = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(4).warped             = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(5).warped             = [0 0];
    matlabbatch{4}.spm.spatial.preproc.tissue(6).warped             = [0 0];
    matlabbatch{4}.spm.spatial.preproc.warp.mrf                     = 1;
    matlabbatch{4}.spm.spatial.preproc.warp.cleanup                 = 1;
    matlabbatch{4}.spm.spatial.preproc.warp.reg                     = [0 0.001 0.5 0.5 0.2];
    matlabbatch{4}.spm.spatial.preproc.warp.affreg                  = 'mni';
    matlabbatch{4}.spm.spatial.preproc.warp.fwhm                    = 0;
    matlabbatch{4}.spm.spatial.preproc.warp.samp                    = [1 1];
    matlabbatch{4}.spm.spatial.preproc.warp.write                   = [1 1];
    matlabbatch{4}.spm.spatial.preproc.channel.vols                 = {[subj_dir,'anat\',subj_handle,'_T1w.nii,1']};
    
    % normalise structural
    matlabbatch{5}.spm.spatial.normalise.write.woptions.bb          = [-78 -112 -70; 78 76 85];
    matlabbatch{5}.spm.spatial.normalise.write.woptions.vox         = [1 1 1];
    matlabbatch{5}.spm.spatial.normalise.write.woptions.interp      = 4;
    matlabbatch{5}.spm.spatial.normalise.write.subj.def             = {[subj_dir,'anat\y_',subj_handle,'_T1w.nii']};
    matlabbatch{5}.spm.spatial.normalise.write.subj.resample        = {[subj_dir,'anat\',subj_handle,'_T1w.nii,1']};
    
    % normalise functional
    matlabbatch{6}.spm.spatial.normalise.write.woptions.bb          = [-78 -112 -70; 78 76 85];
    matlabbatch{6}.spm.spatial.normalise.write.woptions.vox         = [3 3 4];
    matlabbatch{6}.spm.spatial.normalise.write.woptions.interp      = 4;
    matlabbatch{6}.spm.spatial.normalise.write.subj.def             = {[subj_dir,'anat\y_',subj_handle,'_T1w.nii']};
    matlabbatch{6}.spm.spatial.normalise.write.subj.resample        = update_functional_files(matlabbatch{1}.spm.temporal.st.scans,'ua');
    
    % smooth
    matlabbatch{7}.spm.spatial.smooth.fwhm                          = [8 8 8];
    matlabbatch{7}.spm.spatial.smooth.dtype                         = 0;
    matlabbatch{7}.spm.spatial.smooth.im                            = 0;
    matlabbatch{7}.spm.spatial.smooth.prefix                        = 's';
    matlabbatch{7}.spm.spatial.smooth.data                          = update_functional_files(matlabbatch{1}.spm.temporal.st.scans,'wua');
    
    % run batch
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
    % make derivatives folder for subject
    mkdir([data_dir,'derivatives\',subj_handle,'\func\'])
    mkdir([data_dir,'derivatives\',subj_handle,'\anat\'])

    % copy all files to derivatives folder
    copyfile([subj_dir,'func\'],[data_dir,'derivatives\',subj_handle,'\func\'])
    copyfile([subj_dir,'anat\'],[data_dir,'derivatives\',subj_handle,'\anat\'])

    % get functional bids files in derivatives folder
    bids_files = dir([data_dir,'derivatives\',subj_handle,'\func\sub*']);

    % cycle through each file
    for i = 1 : numel(bids_files)

        % delete bids files from derivaties folder
        delete([bids_files(i).folder,'\',bids_files(i).name])
    end
    
    % get anatomical bids files in derivatives folder
    bids_files = dir([data_dir,'derivatives\',subj_handle,'\anat\sub*']);

    % cycle through each file
    for i = 1 : numel(bids_files)

        % delete bids files from derivaties folder
        delete([bids_files(i).folder,'\',bids_files(i).name])
    end
end

% move derivatives folder to cloud
copyfile([data_dir,'derivatives\'],[bear_dir,'derivatives\'])

% delete local copies
rmdir(data_dir,'s')


