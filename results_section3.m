%% RSA
% initialisation
clearvars
clc

% initalise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% define root directory
if ispc;    dir_root = 'Y:/projects/reinstatement_fidelity/';
            dir_tool = 'Y:/projects/general/';
            dir_repos = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
            copylocal = true;
else;       dir_root = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/';
            dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';
            copylocal = false;
end

% add subfunctions
addpath([dir_root,'scripts/subfunctions'])

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
% dilate section 1 mask to fit searchlight size
dilate_mask([dir_root,'bids_data/derivatives/group/rsa-ers/grand_cluster.nii'],scan_search,scan_vox)

% cycle through each subject
for subj = 1 : n_subj
     
    % define subject name
    subj_handle = sprintf('sub-%02.0f',subj);    
    
    % prepare deformation batch
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def        = {[dir_root,'bids_data/derivatives/',subj_handle,'/anat/iy_',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.util.defs.comp{1}.inv.space              = {[dir_root,'bids_data/',subj_handle,'/anat/',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames             = {[dir_root,'bids_data/derivatives/group/rsa-ers/grand_cluster_dilated.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.weight             = {''};
    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr    = {[dir_root,'bids_data/derivatives/',subj_handle,'/masks/']};
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.file           = {[dir_root,'bids_data/derivatives/',subj_handle,'/func/meanua',subj_handle,'_task-rf_run-1_bold.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.preserve           = 0;
    matlabbatch{1}.spm.util.defs.out{1}.push.fwhm               = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.push.prefix             = '';
    
    % run
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
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
    nii_2 = load_untouch_nii([dir_root,'bids_data/derivatives/',subj_handle,'/masks/wgrand_cluster_dilated.nii']);
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
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-rf_rsa-fullVolume.mat'],'scanVec')
    
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
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-rf_rsa-maskedVolume.mat'],'patterns')
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-rf_rsa-mask.mat'],'mask_idx')
    
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
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-rf_rsa-maskedVolume.mat'])
    
    % get mean pattern across all trials and replicate matrix
    meanPattern = repmat(mean(patterns,1),[size(patterns,1), 1]);

    % subtract mean pattern from data
    patterns = patterns - meanPattern;
    
    % save patterns
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-rf_rsa-maskedDemeanedVolume.mat'],'patterns')
    
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
    load([dir_subj,'/rsa-correlation/',subj_handle,'_task-rf_rsa-maskedDemeanedVolume.mat'])
    
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

    % convert stimulus string into value
    for j = 1 : size(stim,1)
        stim{j,4} = find(ismember(conditions,stim{j,2})); %#ok<SAGROW>
    end

    % extract stimulus values and split into A and B
    sv = stim_details.stimulus(stim_details.modality==1);
    Y.sa = sv(1:numel(sv)/2);
    Y.sb = sv(numel(sv)/2+1:end);

    % calculate trial-wise linear discriminant T
    RDM.ldtB = rsa.stat.fisherDiscrTRDM_trainAtestB(X.aa,Y.a,X.bt,Y.b,size(X.aa,2)-8,Y.sb);
    RDM.ldtA = rsa.stat.fisherDiscrTRDM_trainAtestB(X.ba,Y.b,X.at,Y.a,size(X.ba,2)-8,Y.sa);

    % store stimulus values
    RDM.sb = Y.sb;
    RDM.sa = Y.sa;

    % save subject RDM    
    save([dir_subj,'/rsa-correlation/',subj_handle,'_task-rf_rsa-rdm.mat'],'RDM')

    % get mean T
    sRDMs(subj,i,:,:) = mean(cat(3,RDM.ldtB_avg,RDM.ldtA_avg),3); %#ok<SAGROW>
        
    clear Xa Xb Xa_avg Xb_avg Ya Yb
    
    fprintf('Subject %02.0f of %02.0f complete...\n',subj,n_subj)
end

% save average rdms
save([dir_root,'data/combined/rsa/stats/avg_rdms.mat'],'sRDMs')
if ispc; save([dir_repos,'data/searchlight_rdm.mat'],'sRDMs'); end

%% Extract Vector of Similarity Indices for Each Trial
% predefine vector for RSA data
rsa_vec = zeros(n_subj,n_trials);   

% define visual and audio labels
visualLabels = {'eBIKE','eFARM','eUNDERWATER','eWATERMILL'};
audioLabels  = {'eACCORDIAN','eGUITAR','ePIANO','eTRUMPET'};

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subjHandle = sprintf('subj%02.0f',subj);
      
    % load stimulus details
    load([dir_root,'data/fmri/rsa/data/',subjHandle,'/stim.mat'])

    % load RDM    
    load([dir_root,'data/fmri/rsa/data/',subjHandle,'/cluster_221018_RDM.mat'])

    % create cell to hold RSA details
    tmp = cell(n_trials,1);

    % predefine count
    count = 1;
    
    % extract encoding-retrieval pairs
    for j = 1 : size(stim,1)

        % find last occurence of word in 'stim' matrix
        idx = find(ismember(stim(:,1),stim{j,1}),1,'last');

        % if a link exists (only possible if word at encoding)
        if j ~= idx

            % find all matching dynamic stimuli
            match_idx = find(ismember(stim(:,2),stim{j,2}));

            % if matching stimuli are videos
            if any(ismember(visualLabels,stim{j,2}))

                % define mismatching labels
                mm_label = visualLabels(~ismember(visualLabels,stim{j,2}));

                % find mismatching visual stimuli
                mismatch_idx = find(ismember(stim(:,2),mm_label));

            else
                % define mismatching labels
                mm_label = audioLabels(~ismember(audioLabels,stim{j,2}));

                % find mismatching visual stimuli
                mismatch_idx = find(ismember(stim(:,2),mm_label));
            end

            % define which fold data belong to
            if idx <= 192

                % select labels from same fold
                match_idx = match_idx(match_idx <= 192);
                mismatch_idx = mismatch_idx(mismatch_idx <= 192);

                % add rsa value to rsa vector
                tmp{count,1} = mean(RDM.ldtA(idx,match_idx)) - mean(RDM.ldtA(idx,mismatch_idx));%
                
            else         
                % select labels from same fold
                match_idx = match_idx(match_idx > 192) - 192;
                mismatch_idx = mismatch_idx(mismatch_idx > 192) - 192;

                % add rsa value to rsa vector
                tmp{count,1} = mean(RDM.ldtB(idx-192,match_idx)) - mean(RDM.ldtB(idx-192,mismatch_idx));%
            end

            % update count
            count = count + 1;
        end
    end

    % flip rsa values from discriminability to similarity
    rsa_vec(subj,:) = -cell2mat(tmp);
    
    % clear details
    clear count tmp idx j match_idx mm_label mismatch_idx stim RDM subjHandle
end

% save RSA vector
save([dir_root,'data/fmri/rsa/data/similarity_vector.mat'],'rsa_vec')

% clear details
clear visualLabels audioLabels rsa_vec

%% Create AAL Source Grids
% load template grid
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d4mm.mat']); 
template_grid = sourcemodel; clear sourcemodel % rename grid

% load AAL atlas
mri = ft_read_mri([dir_tool,'fieldtrip-20170319/template/atlas/aal/ROI_MNI_V4.nii']);
 
% change MRI to binary 'in-atlas' vs. 'out-atlas'
mri.anatomy = double(mri.anatomy > 0); 

% interpolate mri with grid
cfg                 = [];
cfg.parameter       = 'anatomy';
cfg.interpmethod    = 'nearest';
roi                 = ft_sourceinterpolate(cfg,mri,template_grid);

% clean up
clear mri template_grid cfg

%% Get Time-Frequency Data
% cycle through each subject
for subj = 1 : n_subj
    
    % load
    load([dir_root,'data/eeg/source/subj',sprintf('%02.0f',subj),'_source.mat'])
    
    % get splits
    [audio_bool,hit_bool] = get_splits_eeg(source);
  
    % get time-frequency representation of data
    cfg             = [];
    cfg.trials      = audio_bool==0;
    cfg.keeptrials  = 'yes';
    cfg.method      = 'wavelet';
    cfg.width       = 6;
    cfg.output      = 'pow';
    cfg.toi         = 0.5:0.25:1.5;
    cfg.foi         = 8:2:24; % 50hz sampling
    cfg.pad         = 'nextpow2';
    freq            = ft_freqanalysis(cfg, source); clear source

    % convert powspctrm to single
    freq.powspctrm = single(freq.powspctrm);
    
    % z-transform
    avg_pow = repmat(nanmean(nanmean(freq.powspctrm,4),1),[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
    std_pow = repmat(nanstd(nanmean(freq.powspctrm,4),[],1),[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
    freq.powspctrm = (freq.powspctrm - avg_pow) ./ std_pow; clear avg_pow std_pow
    
    % avg over time and frequency
    cfg             = [];
    cfg.avgovertime = 'yes';
    cfg.avgoverfreq = 'yes';
    freq            = ft_selectdata(cfg,freq);
    
    % save
    save([dir_root,'data/combined/rsa/data/subj',sprintf('%02.0f',subj),'/freq_source.mat'],'freq','-v7.3')
    clear cfg freq    
end

%% Correlate RSA with Source Data
% load RSA data
load([dir_root,'data/fmri/rsa/data/similarity_vector.mat']);

% predefine matrices for chance data
obs_r       = nan(n_subj,2,23156);
chance_r    = nan(n_subj,2,23156);
n_elements  = nan(n_subj,2);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subjHandle = sprintf('subj%02.0f',subj);
      
    % load source data
    load([dir_root,'data/combined/rsa/data/',subjHandle,'/freq_source.mat'])

    % get trl num for eeg data
    trl_num = zeros(size(freq.trialinfo));
    for j = 1 : numel(freq.trialinfo); trl_num(j) = freq.trialinfo{j}.trlNo_encoding; end
    [~,sorted_idx] = sort(trl_num);
        
    % sort powspctrm to match encoding order
    freq.powspctrm = freq.powspctrm(sorted_idx,:);
    freq.trialinfo = freq.trialinfo(sorted_idx,1);
     
    % get splits (in sorted order)
    [~,hit_bool] = get_splits_eeg(freq);

    for i = 1:2
        
        % extract visual trials from RSA matrix
        Y = rsa_vec(subj,ismember(1:192,trl_num))';

        % restrict EEG and RSA data to hits
        X = freq.powspctrm(hit_bool==i-1,:);
        Y = Y(hit_bool==i-1,:);

        % replicate matrix of Y to aid correlation
        Y = repmat(Y,[1,size(X,2)]);

        % correlate
        numer = sum((X-mean(X,1)).*(Y-mean(Y,1)));
        denom = sqrt(sum((X-mean(X,1)).^2)).*sqrt(sum((Y-mean(Y,1)).^2));
        obs_r(subj,i,:) = numer ./ denom;

        % get random correlation
        for k = 1 : 1000

            % reorder X
            Xk = X(randperm(size(X,1)),:);

            % correlate
            numer = sum((Xk-mean(Xk,1)).*(Y-mean(Y,1)));
            denom = sqrt(sum((Xk-mean(Xk,1)).^2)).*sqrt(sum((Y-mean(Y,1)).^2));
            tmp(k,:)   = numer ./ denom;
        end

        % get averaged permutated R
        chance_r(subj,i,:) = mean(tmp);   
        
        % record number of elements in correlation
        n_elements(subj,i) = size(X,1);
    end
    
    % clean up
    clear numer denom tmp Xk X Y k hit_bool freq sorted_idx trl_num j subjHandle
end

% save
save([dir_root,'data/combined/rsa/data/group_corr.mat'],'obs_r','chance_r','n_elements')

%% Run Statistics
% load correlation data
load([dir_root,'data/combined/rsa/data/group_corr.mat'])

% load in source template
load([dir_root,'data/eeg/source/grand_freq.mat'])

% load in EEG cluster
load([dir_root,'data/eeg/source/stat.mat'])

% create frequency structure for analysis
data_hits                               = grand_hits.video;
data_misses                             = grand_hits.video;
data_hits.pow(:,data_hits.inside)       = obs_r(:,2,:);
data_misses.pow(:,data_misses.inside)   = obs_r(:,1,:);
data_misses.pow(:,data_misses.inside)   = obs_r(:,1,:);

% create null hypothesis
null_hits                               = data_hits;
null_misses                             = data_misses;
null_hits.pow(:,null_hits.inside)       = chance_r(:,2,:);
null_misses.pow(:,null_misses.inside)   = chance_r(:,1,:);

% set inside to cluster
data_hits.inside    = stat.video.negclusterslabelmat == 1;
data_misses.inside  = stat.video.negclusterslabelmat == 1;
null_hits.inside    = stat.video.negclusterslabelmat == 1;
null_misses.inside  = stat.video.negclusterslabelmat == 1;
clear stat grand_hits grand_misses

% run statisitics
cfg                   = [];
cfg.dim               = data_hits.dim;  % specify dimensions of your source grid
cfg.design            = [1:size(data_hits.pow,1), 1:size(data_hits.pow,1); ones(1,size(data_hits.pow,1)), ones(1,size(data_hits.pow,1))+1];
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
stat.hits             = ft_sourcestatistics(cfg, data_hits, null_hits);
stat.misses           = ft_sourcestatistics(cfg, data_misses, null_misses);

% save stats
save([dir_root,'data/combined/rsa/stats/eeg_rsa_corr.mat'],'stat')

%% Extract Raw Power of Cluster 
% get indices of clusters
clus_idx = stat.hits.negclusterslabelmat==1;

% get mean magnitude of cluster
hits     = nanmean(data_hits.pow(:,clus_idx),2);
misses   = nanmean(data_misses.pow(:,clus_idx),2);

% create table
tbl = table(hits,misses);

% write table
writetable(tbl,[dir_repos,'data/RSA-EEG_correlation.csv'],'Delimiter',',')

%% Create Surface File
% find indices of largest cluster
clus_idx = stat.hits.negclusterslabelmat==1;

% create source data structure
source                  = [];
source.inside           = stat.hits.inside;
source.dim              = stat.hits.dim;
source.pos              = stat.hits.pos*10;
source.unit             = 'mm';
source.pow              = stat.hits.stat;
source.pow              = nan(size(stat.hits.pos,1),1);     
source.pow(clus_idx)	= stat.hits.stat(clus_idx); 

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
cfg.filename      = 'E:\brain.nii';  % enter the desired file name
cfg.filetype      = 'nifti';
cfg.coordsys      = 'spm';
cfg.vmpversion    = 2;
cfg.datatype      = 'float';
ft_volumewrite(cfg,source);      % be sure to use your interpolated source data


