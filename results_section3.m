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
scan_search = 10;                                           % searchlight radius

% define mask names
mask_names = {'percept','ers'};

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
     
    % cycle through each mask
    for mask = 1 : numel(mask_names)

        % define subject name
        subj_handle = sprintf('sub-%02.0f',subj);    

        % prepare deformation batch
        matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def        = {[dir_root,'bids_data/derivatives/',subj_handle,'/anat/iy_',subj_handle,'_T1w.nii']};
        matlabbatch{1}.spm.util.defs.comp{1}.inv.space              = {[dir_root,'bids_data/',subj_handle,'/anat/',subj_handle,'_T1w.nii']};
        matlabbatch{1}.spm.util.defs.out{1}.push.fnames             = {[dir_root,'bids_data/derivatives/group/rsa-',mask_names{mask},'/grand_cluster_dilated.nii']};
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
        movefile([dir_root,'bids_data/derivatives/',subj_handle,'/masks/wgrand_cluster_dilated.nii'],...
                 [dir_root,'bids_data/derivatives/',subj_handle,'/masks/rsa-',mask_names{mask},'.nii'])
    end
end

%% Read Data
% cycle through each subject
for subj = 1 : n_subj
    
    % update command line
    fprintf('\n--- Working on Subject %d ---------\n',subj)
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/',subj_handle,'/'];
    
    % predefine matrix for functional data
    scanVec = zeros((n_volumes-8).*n_runs,prod(scan_fov));
    
    % start scan counter
    scanCount = 1;
    
    % cycle through each run
    for i = 1 : n_runs
        
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
        fprintf('Run %d of %d read in...\n',i,n_runs)
    end
    
    % clear up
    clear i j filename scanCount nii
    
    % save
    fprintf('\nSaving full volume...\n')
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/'])
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-rf_rsa-fullVolume.mat'],'scanVec')
    
    % apply mask to scans
    fprintf('Masking data...\n')
       
    % predefine cell for masked data
    patterns = cell(numel(mask_names),1);
    
    % cycle through each mask
    for mask = 1 : numel(mask_names)
       
        % predefine matrix for mask data
        maskImg = zeros(1, prod(scan_fov));

        % load mask and add to matrix    
        nii_1 = load_untouch_nii([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/mask.nii']);
        nii_2 = load_untouch_nii([dir_root,'bids_data/derivatives/',subj_handle,'/masks/rsa-',mask_names{mask},'.nii']);
        maskImg(1,:) = reshape(nii_1.img==1&nii_2.img==1,1,[]);

        % define patterns for mask
        patterns{mask} = scanVec;

        % mask functional data (set all zero-element voxels in mask to zero)
        patterns{mask}(:,maskImg(1,:)==0) = 0;

        % get a boolean vector of non-zero voxels in patterns matrix
        mask_idx = all(patterns{mask}~=0);

        % remove zero elements from patterns
        patterns{mask}(:,any(patterns{mask} == 0)) = [];
    end

    % save
    fprintf('Saving masked volumes...\n')
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-all_rsa-maskedVolume.mat'],'patterns')
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-all_rsa-mask.mat'],'mask_idx')

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
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-all_rsa-maskedVolume.mat'])
    
    % cycle through each mask
    for mask = 1 : numel(mask_names)

        % get mean pattern across all trials and replicate matrix
        meanPattern = repmat(mean(patterns{mask},1),[size(patterns{mask},1), 1]);

        % subtract mean pattern from data
        patterns{mask} = patterns{mask} - meanPattern;
    end
    
    % save patterns
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-all_rsa-maskedDemeanedVolume.mat'],'patterns')
    
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
    load([dir_subj,'/rsa-correlation/',subj_handle,'_task-all_rsa-maskedDemeanedVolume.mat'])
    
    % load SPM.mat
    load([dir_subj,'/rsa-correlation/SPM.mat'])
    
    % get stimulus tables
    [stim_details,scan_details] = get_stimulus_tables(dir_root,subj_handle);
    
    % cycle through each mask
    for mask = 1 : numel(mask_names)

        % get design matrix (X)
        X     = [];
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
        Y.raw = patterns{mask}(scan_details.modality==1,:);

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
        clear sv i j

        % reorder from trial order to stimulus order
        [Y.sa,Y.sai] = sort(Y.sa);
        [Y.sb,Y.sbi] = sort(Y.sb);
        X.at(:,1:96) = X.at(:,Y.sai);
        X.bt(:,1:96) = X.bt(:,Y.sbi);

        % switch approach based on data type
        if mask == 1
            % calculate trial-wise linear discriminant T
            RDM.ldtB = rsa.stat.fisherDiscrTRDM_trainAtestB(X.at,Y.a,X.bt,Y.b,size(X.at,2)-96,1:numel(Y.sa));
            RDM.ldtA = rsa.stat.fisherDiscrTRDM_trainAtestB(X.bt,Y.b,X.at,Y.a,size(X.at,2)-96,1:numel(Y.sb));
                        
        else                    
            % calculate trial-wise linear discriminant T
            RDM.ldtB = rsa.stat.fisherDiscrTRDM_trainAtestB(X.aa,Y.a,X.bt,Y.b,size(X.aa,2)-8,Y.sb);
            RDM.ldtA = rsa.stat.fisherDiscrTRDM_trainAtestB(X.ba,Y.b,X.at,Y.a,size(X.ba,2)-8,Y.sa);
            
        end

        % store stimulus values
        RDM.sb = Y.sb;
        RDM.sa = Y.sa;

        % record original trial numbers
        [~,RDM.ta] = sort(Y.sai);
        [~,RDM.tb] = sort(Y.sbi);

        % save subject RDM    
        save([dir_subj,'/rsa-correlation/',subj_handle,'_task-',mask_names{mask},'_rsa-rdm.mat'],'RDM')
        clear X Y RDM
    end

    fprintf('Subject %02.0f of %02.0f complete...\n',subj,n_subj)
end

%% Extract Vector of Similarity Indices for Each Trial
% predefine memory performance
mem_perf = nan(n_subj,96);

% predefine vector for RSA data
rsa_vec  = nan(n_subj,numel(mask_names),n_trials/2);

% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    
    % cycle through each mask
    for mask = 1 : numel(mask_names)
    
        % load RDM
        load([dir_subj,'/rsa-correlation/',subj_handle,'_task-',mask_names{mask},'_rsa-rdm.mat'])

        % switch based on mask
        switch mask
            case 1
                
                % cycle through every encoding trial in RDM
                for trl = 1 : 48

                    % calculate simliarity index
                    rsa_vec(subj,mask,trl)       = calculate_similarity_index(trl,RDM.ldtA,RDM.sa);
                    rsa_vec(subj,mask,trl+48)    = calculate_similarity_index(trl,RDM.ldtB,RDM.sb);
                end

                % reorder trial order from stimulus order to trial order
                rsa_vec(subj,mask,1:48)     = rsa_vec(subj,mask,RDM.ta(1:48));
                rsa_vec(subj,mask,49:96)    = rsa_vec(subj,mask,RDM.tb(1:48)+48);             
                
            case 2     
                
                % cycle through every trial in RDM
                for trl = 49 : size(RDM.ldtA,1)

                    % calculate simliarity index
                    rsa_vec(subj,mask,trl-48)  = calculate_similarity_index(trl,RDM.ldtA,RDM.sa);
                    rsa_vec(subj,mask,trl)     = calculate_similarity_index(trl,RDM.ldtB,RDM.sb);
                end

                % reorder trial order from stimulus order to trial order
                rsa_vec(subj,mask,1:48)     = rsa_vec(subj,mask,RDM.ta(49:end)-48);
                rsa_vec(subj,mask,49:96)    = rsa_vec(subj,mask,RDM.tb(49:end));
        end
        
        % tidy up
        clear RDM trl
        
    end

    % clear details
    clear subj_handle dir_subj
end

% save RSA vector
save([dir_root,'bids_data/derivatives/group/rsa-correlation/group_task-all_fmri-si'],'rsa_vec')

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

% create masks
for mask = 1 : numel(mask_names)

    % load mask
    mask_roi{mask} = ft_read_mri([dir_root,'bids_data/derivatives/group/rsa-',mask_names{mask},'/grand_cluster_dilated.nii']);

    % interpolate
    mask_roi{mask} = ft_sourceinterpolate(cfg,mask_roi{mask},sourcemodel);

    % add
    mask_roi{mask}.inside	= mask_roi{mask}.anatomy(:) > 0 & roi.inside == 1;
    mask_roi{mask}.pos    	= sourcemodel.pos;
end

%% Get Timelocked Representation of Each Condition
% predefine cell for group data
group_freq   = cell(n_subj,2);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];
    
    % load in raw data
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'])  
    
    % get time-frequency representaiton of data for specified conditions
    group_freq{subj,1} = get_basic_timefrequency(source,'encoding','visual');
    group_freq{subj,2} = get_basic_timefrequency(source,'retrieval','visual');
    
    % update command line
    fprintf('\nSubject %02.0f of %02.0f complete...\n\n',subj,n_subj)
end
 
% save data
save([dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-freq.mat'],'group_freq','-v7.3'); 

%% Correlate Similarity and Power
% load EEG stat
load([dir_bids,'derivatives/group/eeg/group_task-all_eeg-stat.mat'])

% load in similarity index
load([dir_root,'bids_data/derivatives/group/rsa-correlation/group_task-all_fmri-si'])

% predefine matrix for correlation values
r = zeros(n_subj,numel(mask_names)+1,2);

% cycle through each subject
for subj = 1 : n_subj
    
    % cycle through each mask (plus one for forgotten ers)
    for mask = 1 : numel(mask_names)+1
    
        % switch based on mask
        switch mask
            case 1
                % get channels in roi
                coi = stat{mask}.negclusterslabelmat(stat{mask}.inside==1)==1;
        
                % extract all trial numbers from powspctrm
                idx         = 1 : size(group_freq{subj,mask}.trialinfo,1);
                trl_nums    = group_freq{subj,mask}.trialinfo(idx,1);
                
            case 2
                % get channels in roi
                coi = stat{mask}.negclusterslabelmat(stat{mask}.inside==1)==1;
                
                % extract recalled trial numbers from powspctrm
                idx         = group_freq{subj,mask}.trialinfo(:,2)==1;
                trl_nums    = group_freq{subj,mask}.trialinfo(idx,1);

            case 3
                % get channels in roi
                coi = stat{mask-1}.negclusterslabelmat(stat{mask-1}.inside==1)==1;
                
                % extract forgotten trial numbers from powspctrm
                idx         = group_freq{subj,mask-1}.trialinfo(:,2)==0;
                trl_nums    = group_freq{subj,mask-1}.trialinfo(idx,1);   
        end

        % fix numbers
        trl_nums(trl_nums>48 & trl_nums<=96)    = trl_nums(trl_nums>48 & trl_nums<=96) - 48;
        trl_nums(trl_nums>96 & trl_nums<=144)   = trl_nums(trl_nums>96 & trl_nums<=144) - 48;
        trl_nums(trl_nums>144 & trl_nums<=192)  = trl_nums(trl_nums>144 & trl_nums<=192) - 96;

        % switch based on mask
        if mask ~= 3

            % extract vectors for items
            X = squeeze(rsa_vec(subj,mask,trl_nums));
            Y = nanmean(nanmean(group_freq{subj,mask}.powspctrm(idx,coi,:),3),2);
            Z = group_freq{subj,mask}.trialinfo(idx,3);
            
        else
            % extract vectors for items
            X = squeeze(rsa_vec(subj,mask-1,trl_nums));
            Y = nanmean(nanmean(group_freq{subj,mask-1}.powspctrm(idx,coi,:),3),2);
            Z = group_freq{subj,mask-1}.trialinfo(idx,3);
        end
        
        % correlate
        r(subj,mask,1) = atanh(corr(X,Y));
        r(subj,mask,2) = atanh(partialcorr(X,Y,Z));
    end
end

% prepare data structure for one-sample t-test
grand_freq{1}              = [];
grand_freq{1}.time         = 1;
grand_freq{1}.freq         = 1;
grand_freq{1}.label        = {'dummy'};
grand_freq{1}.dimord       = 'subj_chan_freq_time';
grand_freq{1}.powspctrm    = r(:,1,1);

% duplicate data structure for ERS recalled data
grand_freq{2}               = grand_freq{1};
grand_freq{2}.powspctrm     = r(:,2,1);

% duplicate data structure for ERS forgotten data
grand_freq{3}               = grand_freq{1};
grand_freq{3}.powspctrm     = r(:,3,1);

% save data
save([dir_bids,'derivatives/group/rsa-correlation/group_task-all_comb-freq.mat'],'grand_freq'); 

%% Run Statistics
% predefine cell for statistics
cfg             = [];
cfg.tail        = -1;
[stat,tbl]      = run_oneSampleT(cfg, grand_freq);
   
% save data
save([dir_bids,'derivatives/group/rsa-correlation/group_task-retrieval_comb-stat.mat'],'stat','tbl');

% check partial correlation controlling for memory confidence
grand_freq_partial = grand_freq{2};
grand_freq_partial.powspctrm = r(:,2,2);

% predefine cell for statistics
cfg             = [];
cfg.tail        = -1;
[stat,tbl]      = run_oneSampleT(cfg, grand_freq_partial);
   
% save data
save([dir_bids,'derivatives/group/rsa-correlation/group_task-retrieval_comb-partialstat.mat'],'stat','tbl');
   
%% Extract Raw Power of Cluster 
% prepare table for stat values
tbl = array2table(zeros(n_subj,4),'VariableNames',{'perception','retrieval','forgotten','partial'});

% create table
tbl.perception(:,1) = grand_freq{1}.powspctrm;
tbl.retrieval(:,1)  = grand_freq{2}.powspctrm;
tbl.forgotten(:,1)  = grand_freq{3}.powspctrm;
tbl.partial(:,1)    = grand_freq_partial.powspctrm;

% write table
writetable(tbl,[dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-cluster.csv'],'Delimiter',',')
copyfile([dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-cluster.csv'],[dir_repos,'data/fig3_data/group_task-all_eeg-cluster.csv'])

%% Get Time/Freq Series of Data
% load EEG stat
load([dir_bids,'derivatives/group/eeg/group_task-all_eeg-stat.mat'])

% load in similarity index
load([dir_root,'bids_data/derivatives/group/rsa-correlation/group_task-all_fmri-si'])

% predefine cell for group data
group_freq   = nan(n_subj,75,81,2);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];
    
    % load in raw data
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'])  
    
    % get channels in roi
    coi = stat{2}.negclusterslabelmat(stat{2}.inside==1)==1;

    % select region of interest
    cfg = [];
    cfg.channel = source.label(coi);
    data = ft_selectdata(cfg,source); clear source
        
    % predefine conditional arrays to include all trials
    operation_to_include = zeros(numel(data.trial),1);
    modality_to_include  = zeros(numel(data.trial),1);

    % predefine arrays for memory performance
    mem_performance = zeros(numel(data.trial),1);
    trl_nums = zeros(numel(data.trial),1);

    % cycle through each trial
    for trl = 1 : numel(data.trial)
        
        % mark trials that do match specified operation
        operation_to_include(trl) = strcmpi(data.trialinfo{trl}.operation,'retrieval');
        
        % mark trials that do match specified operation
        modality_to_include(trl) = strcmpi(data.trialinfo{trl}.modality,'visual');
        
        % get memory
        mem_performance(trl) = data.trialinfo{trl}.recalled;
        trl_nums(trl) = data.trialinfo{trl}.trl_at_retrieval;
    end

    % select data
    cfg                 = [];
    cfg.trials          = operation_to_include == 1 & modality_to_include == 1;
    data                = ft_selectdata(cfg,data);

    % select memory performance details
    mem_performance = mem_performance(operation_to_include == 1 & modality_to_include == 1);
    trl_nums = trl_nums(operation_to_include == 1 & modality_to_include == 1);

    % get time-frequency parameters
    cfg                 = [];
    cfg.keeptrials      = 'yes';
    cfg.method          = 'wavelet';
    cfg.width           = 6;
    cfg.output          = 'pow';
    cfg.toi             = -1:0.05:3;
    cfg.foi             = 3:0.5:40;
    cfg.pad             = 'nextpow2';
    freq                = ft_freqanalysis(cfg, data);
    clear data

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
    
    for i = 1 : 2

        % extract recalled trial numbers from powspctrm
        idx     = mem_performance==i-1;
        trli    = trl_nums(idx,1);

        % fix numbers
        trli(trli>48 & trli<=96)    = trli(trli>48 & trli<=96) - 48;
        trli(trli>96 & trli<=144)   = trli(trli>96 & trli<=144) - 48;
        trli(trli>144 & trli<=192)  = trli(trli>144 & trli<=192) - 96;

        % extract vectors for items
        X = squeeze(rsa_vec(subj,2,trli));
        
        % cycle through time nad freq
        for t = 1 : size(freq.powspctrm,4)
            for f = 1 : size(freq.powspctrm,3)
        
                Y = nanmean(freq.powspctrm(idx,:,f,t),2);

                group_freq(subj,f,t,i) = corr(X,Y);
            end
        end
    end

    % update command line
    fprintf('\nSubject %02.0f of %02.0f complete...\n\n',subj,n_subj)
end
 
% define time and frequency vectors
t = -1:0.05:3;
f = 3:0.5:40;

% extract time and frequency matrices
raw_time = squeeze(mean(group_freq(:,f>=8&f<=30,:,:),2));
raw_freq = squeeze(mean(group_freq(:,:,t>=0.5&t<=1.5,:),3));

% create table for time data
tbl_time = array2table([raw_time(:),...
                        repmat((1:n_subj)',[numel(raw_time)/n_subj,1]),...
                        cat(1,reshape(repmat(t,[numel(raw_time)/numel(t)/2,1]),[],1),reshape(repmat(t,[numel(raw_time)/numel(t)/2,1]),[],1)),...
                        [ones(numel(raw_time)/2,1);ones(numel(raw_time)/2,1)+1]],...
                        'VariableNames',{'signal','subject','time','condition'});

% create table for freq data
tbl_freq = array2table([raw_freq(:),...
                        repmat((1:n_subj)',[numel(raw_freq)/n_subj,1]),...
                        cat(1,reshape(repmat(f,[numel(raw_freq)/numel(f)/2,1]),[],1),reshape(repmat(f,[numel(raw_freq)/numel(f)/2,1]),[],1)),...
                        [ones(numel(raw_freq)/2,1);ones(numel(raw_freq)/2,1)+1]],...
                        'VariableNames',{'signal','subject','freq','condition'});
          
% write table
writetable(tbl_time,[dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-timeseries.csv'],'Delimiter',',')
writetable(tbl_freq,[dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-freqseries.csv'],'Delimiter',',')

% copy to repository
copyfile([dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-timeseries.csv'],[dir_repos,'data/fig3_data/group_task-all_eeg-timeseries.csv'])
copyfile([dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-freqseries.csv'],[dir_repos,'data/fig3_data/group_task-all_eeg-freqseries.csv'])

