%% RSA
% initialisation
clearvars
clc

% initalise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% define root directory
if ispc;    dir_root = 'Y:/projects/reinstatement_fidelity/';
            dir_bids = 'Y:/projects/reinstatement_fidelity/bids_data/';
            dir_tool = 'Y:/projects/general/';
            dir_repos = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
            copylocal = true;
else;       dir_root = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/';
            dir_bids = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/bids_data/';
            dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';
            copylocal = false;
end

% add subfunctions
addpath([dir_root,'scripts/subfunctions'])

% add developer rsa toolbox
addpath([dir_tool,'rsatoolbox-develop'])

% load layout
load([dir_root,'data/layout.mat'])

%% Define Key Parameters
dir_root    = 'Y:/projects/reinstatement_fidelity/';        % data directory
dir_repos   = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
n_subj      = 21;                                           % number of subjects
n_trials    = 192;                                          % number of trials
n_volumes   = 255;                                          % number of volumes
n_runs      = 8;
TR          = 2;
EEG_sample  = 5000;
scan_fov    = [64 64 32];                                   % scan field of view
scan_vox    = [3 3 4];                                      % scan voxel size
scan_search = 8;                                           % searchlight radius
scan_func   = {'_3_1','_4_1','_5_1','_6_1',...
               '_8_1','_9_1','_10_1','_11_1'};              % functional scan suffix

% add subfunctions
addpath([dir_repos,'subfunctions'])

%% Create GLM Model
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    mkdir([dir_subj,'rsa-correlation/'])
    
    % delete old SPM file if it exists
    if exist([dir_subj,'rsa-correlation/SPM.mat'],'file'); delete([dir_subj,'rsa-correlation/SPM.mat']); end
    
    % create cell to hold events data
    events_onset  = zeros(n_trials/2,1);
    count = 1;
    
    % create array to hold button press onsets
    button_onset = [];
    
    % define 'length of scan'
    LoS = 0;
        
    % cycle through each run
    for run = 1 : n_runs

        % load event table
        tbl = readtable([dir_root,'bids_data/',subj_handle,'/func/',subj_handle,'_task-rf_run-',num2str(run),'_events.tsv'],'FileType','text','Delimiter','\t');
        
        % remove first three volumes and last five volumes
        tbl = tbl(4:end-5,:);
        
        % adjust time accordingly
        tbl.onset = tbl.onset - 30000;
        
        % add length of prior scans to this event table
        tbl.onset = tbl.onset + LoS;
        
        % extract 'length of scan' and add on additional scan for next iteration
        LoS = max(tbl.onset) + (TR*EEG_sample)-1;
        
        % cut table to stimulus onset
        tbl = tbl(ismember(tbl.trial_type,'Stimulus Onset'),:);
        
        % cycle through every event
        for e = 1 : size(tbl,1)
            
            % check if encoding adn visual
            if strcmpi(tbl.modality{e},'visual')

                % add key values
                events_onset(count,1) = tbl.onset(e);
                count = count + 1;
            end

            % get button press onset (if pressed)
            if ~isnan(tbl.rt(e)); button_onset(end+1,1) = tbl.onset(e) + tbl.rt(e); end %#ok<SAGROW>            
        end
    end
    
    % tidy up
    clear count e tbl LoS run
       
    % convert event onsets and button presses from EEG samples to seconds
    events_onset = events_onset ./ EEG_sample;
    button_onset = button_onset ./ EEG_sample;    
    
    % load movement parameters
    R = load([dir_subj,'/func/rp_a',subj_handle,'_task-rf_run-1_bold.txt']);
    
    % find first three volumes and last five volumes of each run
    bad_scans = zeros(8*n_runs,1);
    for run = 1 : n_runs
        bad_scans(1+(8*(run-1)):8*run) = [1 2 3 251 252 253 254 255] + n_volumes*(run-1);
    end
    
    % remove these scans
    R(bad_scans,:) = [];
    
    % define number of volumes remaining in a run
    n_volumes_adj = n_volumes - 8;
    
    % add regressors that model the per-block linear change
    R(1+(n_volumes_adj*0):n_volumes_adj*1,7) = linspace(0,1,n_volumes_adj);
    R(1+(n_volumes_adj*2):n_volumes_adj*3,8) = linspace(0,1,n_volumes_adj);
    R(1+(n_volumes_adj*4):n_volumes_adj*5,9) = linspace(0,1,n_volumes_adj);
    R(1+(n_volumes_adj*6):n_volumes_adj*7,10) = linspace(0,1,n_volumes_adj);
    R(1+(n_volumes_adj*1):n_volumes_adj*2,11) = linspace(0,1,n_volumes_adj);
    R(1+(n_volumes_adj*3):n_volumes_adj*4,12) = linspace(0,1,n_volumes_adj);
    R(1+(n_volumes_adj*5):n_volumes_adj*6,13) = linspace(0,1,n_volumes_adj);
    R(1+(n_volumes_adj*7):n_volumes_adj*8,14) = linspace(0,1,n_volumes_adj);
    
    % add regressors that model the per-block constant
    R(1+(n_volumes_adj*0):n_volumes_adj*1,15) = ones(1,n_volumes_adj);
    R(1+(n_volumes_adj*2):n_volumes_adj*3,16) = ones(1,n_volumes_adj);
    R(1+(n_volumes_adj*4):n_volumes_adj*5,17) = ones(1,n_volumes_adj);
    R(1+(n_volumes_adj*6):n_volumes_adj*7,18) = ones(1,n_volumes_adj);
    R(1+(n_volumes_adj*1):n_volumes_adj*2,19) = ones(1,n_volumes_adj);
    R(1+(n_volumes_adj*3):n_volumes_adj*4,20) = ones(1,n_volumes_adj);
    R(1+(n_volumes_adj*5):n_volumes_adj*6,21) = ones(1,n_volumes_adj);
    R(1+(n_volumes_adj*7):n_volumes_adj*8,21) = ones(1,n_volumes_adj); 
    
    % save nuisance regressors
    if copylocal; mkdir('C:/tmp_mri/spm/'); save('C:/tmp_mri/spm/R.mat','R')   
    else; save([dir_subj,'/rsa-correlation/R.mat'],'R')   
    end
              
    clear R n_volumes_adj run
    
    % get all scans for GLM
    all_scans = get_functional_files(dir_subj,'ua',copylocal);
    all_scans = all_scans{1};
    
    % remove bad scans
    all_scans(bad_scans) = [];
    
    % define parameters for GLM
    matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
    matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi       = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress     = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf         = 128;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans       = all_scans;  
    matlabbatch{1}.spm.stats.fmri_spec.timing.units     = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT        = 2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t    = 32;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0   = 16;  
    
    % define directories
    if copylocal
        matlabbatch{1}.spm.stats.fmri_spec.dir              = {'C:/tmp_mri/spm'};
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {'C:/tmp_mri/spm/R.mat'};         
    else
        matlabbatch{1}.spm.stats.fmri_spec.dir              = {[dir_subj,'/rsa-correlation']};
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {[dir_subj,'/rsa-correlation/R.mat']}; 
    end
                   
    % cycle through and define each condition
    for trl = 1 : size(events_onset,1)
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(trl).name        = ['trl',sprintf('%03.0f',trl)];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(trl).onset       = events_onset(trl,1);
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(trl).duration    = 3;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(trl).tmod        = 0;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(trl).orth        = 1;
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(trl).pmod        = struct('name',{},'param',{},'poly',{});
    end
    
    % add button press
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(end+1).name      = 'button_press';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(end).onset       = button_onset;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(end).duration    = 0.5;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(end).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(end).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(end).pmod        = struct('name',{},'param',{},'poly',{});
        
    % estimate model
    matlabbatch{2}.spm.stats.fmri_est.write_residuals                   = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical                  = 1;

    % define directories
    if copylocal; matlabbatch{2}.spm.stats.fmri_est.spmmat = {'C:/tmp_mri/spm/SPM.mat'};
    else; matlabbatch{2}.spm.stats.fmri_est.spmmat = {[dir_subj,'/rsa-correlation/SPM.mat']};
    end
    
    % run batch
    spm_jobman('run',matlabbatch)
    clear matlabbatch all_scans bad_scans subj_handle button_onset events_onset trl    
    
    % copy and delete local files if they exist
    if copylocal
        copyfile('C:/tmp_mri/spm/R.mat',[dir_subj,'/rsa-correlation/R.mat'])
        copyfile('C:/tmp_mri/spm/SPM.mat',[dir_subj,'/rsa-correlation/SPM.mat'])
        copyfile('C:/tmp_mri/spm/mask.nii',[dir_subj,'/rsa-correlation/mask.nii'])
        rmdir('C:/tmp_mri/','s'); 
    end
end

%% Prepare Masks
% cycle through each subject
for subj = 1 : n_subj
     
    % define subject name
    subj_handle = sprintf('sub-%02.0f',subj);    
    
    % prepare deformation batch
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def        = {[dir_root,'bids_data/derivatives/',subj_handle,'/anat/iy_',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.util.defs.comp{1}.inv.space              = {[dir_root,'bids_data/',subj_handle,'/anat/',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames             = {[dir_root,'bids_data/derivatives/group/rsa-percept/grand_cluster.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.weight             = {''};
    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr    = {[dir_root,'bids_data/derivatives/',subj_handle,'/masks/']};
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.file           = {[dir_root,'bids_data/derivatives/',subj_handle,'/func/meanua',subj_handle,'_task-rf_run-1_bold.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.preserve           = 0;
    matlabbatch{1}.spm.util.defs.out{1}.push.fwhm               = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.push.prefix             = '';
    
    % run
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
    % change mask name
    movefile([dir_root,'bids_data/derivatives/',subj_handle,'/masks/wgrand_cluster.nii'],...
             [dir_root,'bids_data/derivatives/',subj_handle,'/masks/rsa-percept.nii'])
    
end

%% Read Data
% cycle through each subject
for subj = 1 : n_subj
    
    % update command line
    fprintf('\n--- Working on Subject %d ---------\n',subj)
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/',subj_handle,'/'];
    
    % predefine matrix for mask data
    maskImg = zeros(1, prod(scan_fov));
    
    % load mask and add to matrix    
    nii_1 = load_untouch_nii([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/mask.nii']);
    nii_2 = load_untouch_nii([dir_root,'bids_data/derivatives/',subj_handle,'/masks/rsa-percept.nii']);
    maskImg(1,:) = reshape(nii_1.img==1&nii_2.img==1,1,[]);

    % predefine matrix for functional data
    scanVec = zeros((n_volumes-8).*numel(scan_func),numel(nii_1.img));
    
    % start scan counter
    scanCount = 1;
    
    % cycle through each run
    for i = 1 : numel(scan_func)
        
        % define functional filenames
        filename = [dir_root,'bids_data/derivatives/',subj_handle,'/func/',...
            'ua',subj_handle,'_task-rf_run-',num2str(i),'_bold.nii'];

        % read in nifti file
        nii = load_untouch_nii(filename);
        
        % define good scans
        good_scans = 4 : n_volumes-5;
        
        % cycle through each scan (excluding first three and last five)
        for j = 1 : numel(good_scans)
            
            % extract image
            scanVec(scanCount,:) = reshape(nii.img(:,:,:,good_scans(j)),1,[]);
            
            % count scan as read in
            scanCount = scanCount + 1;            
        end
        
        % update command line
        fprintf('Run %d of %d read in...\n',i,numel(scan_func))
    end
    
    % clear up
    clear i j filename scanCount nii
    
    % save
    fprintf('\nSaving full volume...\n')
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/'])
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-percept_rsa-fullVolume.mat'],'scanVec')
    
    % apply mask to scans
    fprintf('Masking data...\n')
    
    % define patterns for mask
    patterns = scanVec;

    % mask functional data (set all zero-element voxels in mask to zero)
    patterns(:,maskImg(1,:)==0) = 0;

    % get a boolean vector of non-zero voxels in patterns matrix
    mask_idx = all(patterns~=0);
    
    % remove zero elements from patterns
    patterns(:,any(patterns == 0)) = [];
    
    % save
    fprintf('Saving masked volumes...\n')
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-percept_rsa-maskedVolume.mat'],'patterns')
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-percept_rsa-mask.mat'],'mask_idx')
    
    % clear excess variables
    clear scanVec nii deadIdx patterns mask_idx subjHandle maskImg
end

%% Mean Pattern Subtraction
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/',subj_handle,'/'];
    
    % load pattern data
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-percept_rsa-maskedVolume.mat'])
    
    % get mean pattern across all trials and replicate matrix
    meanPattern = repmat(mean(patterns,1),[size(patterns,1), 1]);

    % subtract mean pattern from data
    patterns = patterns - meanPattern;
    
    % save patterns
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-percept_rsa-maskedDemeanedVolume.mat'],'patterns')
    
    % clean up
    clear subjHandle meanPattern patterns
end

%% Calculate LDt
% cycle through each subject
for subj = 1 : n_subj
       
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    
    % load pattern data
    load([dir_subj,'/rsa-correlation/',subj_handle,'_task-percept_rsa-maskedDemeanedVolume.mat'])
    
    % load SPM.mat
    load([dir_subj,'/rsa-correlation/SPM.mat'])
    
    % get stimulus tables
    [stim_details,scan_details] = get_stimulus_tables(dir_root,subj_handle);
    
    % get design matrix (X)
    X.raw = SPM.xX.X;
        
    % remove scans/regressors that are not visual (to
    % computationally demanding to anything more than this)
    X.raw = X.raw(scan_details.modality==1,:);
    
    % get condition-averaged matrix
    for i = 1 : 8
        X.c(:,i) = sum(X.raw(:,find(stim_details.stimulus(stim_details.modality==1)==i)),2); %#ok<FNDSB>
    end
    
    % add nuisance regressors to condition-averaged matrix
    nR = size(X.raw,2)-n_trials;
    X.c(:,end+1:end+nR) = X.raw(:,end-(nR-1):end);
           
    % split GLM into two groups (train and test data) [excluding nuisance regressors)
    X.at = X.raw(1:size(X.raw,1)/2,[1:n_trials/2 n_trials+1:end]);
    X.bt = X.raw((size(X.raw,1)/2)+1:end,n_trials/2+1:end);
    
    % repeat for condition-average GLM
    X.aa = X.c(1:size(X.raw)/2,:);
    X.ba = X.c((size(X.raw)/2)+1:end,:);
    
    % find zero-value rows
    X.a_badRow = all(X.at(:,1:n_trials/2)==0,2);
    X.b_badRow = all(X.bt(:,1:n_trials/2)==0,2);
    
    % remove zero-value rows
    X.at(X.a_badRow,:) = [];
    X.aa(X.a_badRow,:) = [];
    X.bt(X.b_badRow,:) = [];
    X.ba(X.b_badRow,:) = [];
    
    % cycle through A/B avg/trl combinations
    fn = {'at','aa','bt','ba'};
    for i = 1 : 4
        for j = 1 : size(X.(fn{i}),2)
            
            % if column is singular
            if numel(unique(X.(fn{i})(:,j))) == 1
                X.([fn{i},'_badCol'])(j,1) = true;
            else
                X.([fn{i},'_badCol'])(j,1) = false;
            end
        end
    end
    
    % remove single-valued columns
    X.at(:,X.at_badCol) = [];
    X.aa(:,X.aa_badCol) = [];
    X.bt(:,X.bt_badCol) = [];
    X.ba(:,X.ba_badCol) = [];
      
    % get activation matrix (Y) and 
    Y.raw = patterns(scan_details.modality==1,:);
    
    % split into two groups (Ya and Yb)
    Y.a = Y.raw(1:size(Y.raw,1)/2,:);
    Y.b = Y.raw(size(Y.raw,1)/2+1:end,:);
     
    % remove zero-value rows
    Y.a(X.a_badRow,:) = [];
    Y.b(X.b_badRow,:) = [];

    % extract stimulus values and split into A and B
    sv = stim_details.stimulus(stim_details.modality==1);
    Y.sa = sv(1:numel(sv)/2);
    Y.sb = sv(numel(sv)/2+1:end);
    
    % reorder stimuli
    [Y.sa,Y.sai] = sort(Y.sa);
    [Y.sb,Y.sbi] = sort(Y.sb);
    X.at(:,1:96) = X.at(:,Y.sai);
    X.bt(:,1:96) = X.bt(:,Y.sbi);
    
    % tidy up    
    clear sv i j nR stim_details scan_details SPM patterns

    % calculate trial-wise linear discriminant T
    RDM.ldtB = rsa.stat.fisherDiscrTRDM_trainAtestB(X.aa,Y.a,X.bt,Y.b,size(X.aa,2)-8,Y.sb);
    RDM.ldtA = rsa.stat.fisherDiscrTRDM_trainAtestB(X.ba,Y.b,X.at,Y.a,size(X.ba,2)-8,Y.sa);

    % store stimulus values
    RDM.sb = Y.sb;
    RDM.sa = Y.sa;

    % save subject RDM    
    save([dir_subj,'/rsa-correlation/',subj_handle,'_task-percept_rsa-rdm.mat'],'RDM')
    clear X Y RDM
    
    fprintf('Subject %02.0f of %02.0f complete...\n',subj,n_subj)
end

%% Extract Vector of Similarity Indices for Each Trial
% predefine vector for RSA data
rsa_vec  = zeros(n_subj,n_trials/2);

% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    
    % load RDM
    load([dir_subj,'/rsa-correlation/',subj_handle,'_task-percept_rsa-rdm.mat'])
    
    % cycle through every trial in RDM
    for trl = 1 : size(RDM.ldtA,1)/2
        
        % calculate simliarity index
        rsa_vec(subj,trl)       = calculate_similarity_index(trl,RDM.ldtA,RDM.sa);
        rsa_vec(subj,trl+48)    = calculate_similarity_index(trl,RDM.ldtB,RDM.sb);
    end
    
    % clear details
    clear subj_handle dir_subj RDM trl
end

% save RSA vector
save([dir_root,'bids_data/derivatives/group/rsa-correlation/group_task-percept_fmri-si'],'rsa_vec')

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
group_freq   = cell(n_subj,1);

% load in similarity index
load([dir_root,'bids_data/derivatives/group/rsa-correlation/group_task-percept_fmri-si'])

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];
    
    % load in raw data
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'])  
    
    % get time-frequency representaiton of data for specified conditions
    group_freq{subj,1} = get_rsa_timefrequency_correlation(source,'encoding','visual',rsa_vec(subj,:));
end
   
% get grand average of subjects
cfg                 = [];
cfg.keepindividual  = 'yes';
grand_freq          = ft_freqgrandaverage(cfg,group_freq{:,1});

% collapse over time/frequency
grand_freq.pow = mean(mean(grand_freq.powspctrm,4),3);

% add in source model details
grand_freq.powdimord   = 'rpt_pos';
grand_freq.dim         = roi.dim;
grand_freq.inside      = roi.inside(:);
grand_freq.pos         = roi.pos;
grand_freq.cfg         = [];

% add in "outside" virtual electrods
tmp_pow = zeros(size(grand_freq.pow,1),size(grand_freq.inside,1));
tmp_pow(:,grand_freq.inside) = grand_freq.pow;
grand_freq.pow = tmp_pow;

% remove freq details
grand_freq = rmfield(grand_freq,{'powspctrm','label','freq','time','dimord'});
    
% load EEG stat
load([dir_bids,'derivatives/group/eeg/group_task-all_eeg-stat.mat']);

% switch inside to EEG cluster
grand_freq.inside = stat{1}.negclusterslabelmat == 1;

% save data
save([dir_bids,'derivatives/group/rsa-correlation/group_task-percept_comb-freq.mat'],'grand_freq'); 

%% Run Statistics
% predefine cell for statistics
cfg.tail        = -1;
cfg.parameter   = 'pow';
[stat,tbl]      = run_oneSampleT(cfg, grand_freq);

% save data
save([dir_bids,'derivatives/group/rsa-correlation/group_task-retrieval_comb-stat.mat'],'stat','tbl');


