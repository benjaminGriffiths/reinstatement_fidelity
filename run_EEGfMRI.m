function run_EEGfMRI(skipMRI,skipEEG)

% initalise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% check inputs
if ~exist('skipMRI','var'); skipMRI = true; end
if ~exist('skipEEG','var'); skipEEG = true; end

% define root directory
if ispc;    dir_git  = 'E:/bjg335/projects/reinstatement_fidelity/';
            dir_bids = 'Y:/projects/reinstatement_fidelity/';
            dir_tool = 'Y:/projects/general/';
            dir_repos = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
    
            copylocal = true;
else;       dir_root = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/';
            dir_bids = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/bids_data/';
            dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';
            copylocal = false;
end

% add subfunctions
addpath([dir_git,'/subfunctions'])

% add developer rsa toolbox
addpath([dir_tool,'rsatoolbox-develop'])
addpath([dir_tool,'wen_spectral'])

%% Define Key Parameters
n_subj      = 21;                                           % number of subjects
n_trials    = 192;                                          % number of trials
n_volumes   = 255;                                          % number of volumes
n_runs      = 8;
TR          = 2;
EEG_sample  = 5000;
scan_fov    = [64 64 32];                                   % scan field of view

% define mask names
mask_names = {'percept','ers'};

%% --- RUN FMRI ANALYSIS ----------------------------------------------- %%
if ~skipMRI

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
                if ~isnan(tbl.rt(e)); button_onset(end+1,1) = tbl.onset(e) + tbl.rt(e); end            %#ok<AGROW>
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
        
        % load pattern data
        load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-correlation/',subj_handle,'_task-all_rsa-maskedVolume.mat'],'patterns')

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

    %% Get Trial BOLD 
    % predefine mean bold
    mean_bold = nan(n_subj,numel(mask_names),n_trials);

    % cycle through each subject
    for subj = 1 : n_subj

        % define key subject strings
        subj_handle = sprintf('sub-%02.0f',subj);
        dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];

        % load pattern data
        load([dir_subj,'/rsa-correlation/',subj_handle,'_task-all_rsa-maskedDemeanedVolume.mat'],'patterns')

        % load SPM.mat
        load([dir_subj,'/rsa-correlation/SPM.mat'],'SPM')

        % cycle through each mask
        for mask = 1 : numel(mask_names)

            % get design matrix (X) and patterns (Y)
            X = SPM.xX.X;
            Y = patterns{mask};

            % cycle through each trial
            for trl = 1 : n_trials

                % get scans invovled in trial
                idx = X(:,trl) > 0;
                x = repmat(X(idx,trl),[1 size(Y,2)]);
                y = Y(idx,:);

                % get mean bold signal
                mean_bold(subj,mask,trl) = mean(mean(x .* y,2),1);
            end
        end
    end

    % save
    save([dir_root,'bids_data/derivatives/group/rsa-correlation/group_task-all_fmri-meanbold'],'mean_bold')
    clear X Y x y

    %% Calculate LDt
    % cycle through each subject
    for subj = 1 : n_subj

        % define key subject strings
        subj_handle = sprintf('sub-%02.0f',subj);
        dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];

        % load pattern data
        load([dir_subj,'/rsa-correlation/',subj_handle,'_task-all_rsa-maskedDemeanedVolume.mat'],'patterns')

        % load SPM.mat
        load([dir_subj,'/rsa-correlation/SPM.mat'],'SPM')

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

            % add memory details
            mem = stim_details.memory(stim_details.modality==1);
            Y.ma = mem(1:numel(mem)/2,1);
            Y.mb = mem(numel(mem)/2+1:end,1);

            % reorder from trial order to stimulus order
            [Y.sa,Y.sai] = sort(Y.sa);
            [Y.sb,Y.sbi] = sort(Y.sb);
            X.at(:,1:96) = X.at(:,Y.sai);
            X.bt(:,1:96) = X.bt(:,Y.sbi);
            Y.ma         = Y.ma(Y.sai,:);
            Y.mb         = Y.mb(Y.sbi,:);

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

            % store memory details
            RDM.ma = Y.ma;
            RDM.mb = Y.mb;

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
            load([dir_subj,'/rsa-correlation/',subj_handle,'_task-',mask_names{mask},'_rsa-rdm.mat'],'RDM')

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

                    % extract memory details
                    mem_perf(subj,1:48) = RDM.ma(RDM.ta(49:end));
                    mem_perf(subj,49:96) = RDM.ma(RDM.ta(49:end));
            end

            % tidy up
            clear RDM trl        
        end

        % clear details
        clear subj_handle dir_subj
    end

    mem_rsa = squeeze(rsa_vec(:,2,:));
    figure; hold on
    histogram(mem_rsa(mem_perf==1),-0.975:0.05:0.975)
    histogram(mem_rsa(mem_perf==0),-0.975:0.05:0.975)

    % save RSA vector
    save([dir_root,'bids_data/derivatives/group/rsa-correlation/group_task-all_fmri-si'],'rsa_vec')
    
end

% --- FMRI DONE -------------------------------------------------------- %%

%% Prepare ROI
% load sourcemodel
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat'],'sourcemodel');

% load whole brain
mri = ft_read_mri([dir_bids,'sourcedata/masks/whole_brain.nii']);
 
% load in raw data of subject 1
dir_subj = [dir_bids,'sourcedata/sub-01/eeg/'];
load([dir_subj,'sub-01_task-rf_eeg-source.mat'],'source')  
    
% interpolate clusters with grid
cfg             = [];
cfg.parameter	= 'anatomy';
roi             = ft_sourceinterpolate(cfg,mri,sourcemodel);

% define additional roi parameters
roi.inside      = roi.anatomy(:) > 0;
roi.pos         = sourcemodel.pos;

% create masks
mask_roi = cell(numel(mask_names),1);
for mask = 1 : numel(mask_names)

    % load mask
    mask_roi{mask} = ft_read_mri([dir_bids,'/derivatives/group/rsa-',mask_names{mask},'/grand_cluster_dilated.nii']);

    % interpolate
    mask_roi{mask} = ft_sourceinterpolate(cfg,mask_roi{mask},sourcemodel);
    mask_roi{mask}.anatomy = mask_roi{mask}.anatomy(:);

    % add
    mask_roi{mask}.freq_2_roi = find(mask_roi{mask}.anatomy(roi.inside==1)>0);
    mask_roi{mask}.inside	= mask_roi{mask}.anatomy(:) > 0 & roi.inside == 1;
    mask_roi{mask}.label    = source.label(mask_roi{mask}.inside(roi.inside == 1));
    mask_roi{mask}.pos    	= sourcemodel.pos;
end

% copy complete roi into mask_roi
mask_roi{3} = roi;
mask_roi{4} = mri;

% clean workspace
clear mri source dir_subj roi cfg lay

%% --- RUN EEG ANALYSIS ------------------------------------------------ %%
if ~skipEEG;tic

    %% Get Wavelet Power
    % cycle through each subject
    for subj = 1 : n_subj

        % load in raw data
        fprintf('\nloading sub-%02.0f data...\n',subj);
        load(sprintf('%sderivatives/sub-%02.0f/eeg/sub-%02.0f_task-rf_eeg-wavelet.mat',dir_bids,subj,subj),'freq')  
        
        % recode trialinfo
        freq = recode_trlinfo(freq);
        
        % select encoding data of interest
        cfg             = [];
        cfg.channel     = mask_roi{1}.label;
        cfg.trials      = freq.trialinfo(:,5) == 1 & freq.trialinfo(:,6) == 1;
        cfg.avgovertime = 'yes';
        cfg.avgoverchan = 'yes';
        cfg.latency     = [0.5 1.5];
        tmp{1,1}        = ft_selectdata(cfg,freq);
        cfg.latency     = [-1 -0.375];
        tmp{1,2}        = ft_selectdata(cfg,freq);
        
        % select retrieval data of interest
        cfg.channel     = mask_roi{2}.label;
        cfg.trials      = freq.trialinfo(:,5) == 0 & freq.trialinfo(:,6) == 1;
        cfg.latency     = [0.5 1.5];
        tmp{2,1}        = ft_selectdata(cfg,freq);
        cfg.latency     = [-1 -0.375];
        tmp{2,2}        = ft_selectdata(cfg,freq);
        
        % redefine freq and get ERD
        freq = tmp(:,1);
        freq{1}.powspctrm = cat(4,tmp{1,1}.powspctrm,tmp{1,2}.powspctrm);
        freq{2}.powspctrm = cat(4,tmp{2,1}.powspctrm,tmp{2,2}.powspctrm);
        
        % save spectral outputs
        save(sprintf('%sderivatives/sub-%02.0f/eeg/sub-%02.0f_task-corr_eeg-wavelet.mat',dir_bids,subj,subj),'freq')  
        fprintf('\nsub-%02.0f complete...\ntotal time elapsed: %1.1f minutes...\nestimated time remaining: %1.1f minutes...\n',subj,round(toc/60,1),round((round(toc/60,1)/subj)*(n_subj-subj),1))

        % tidy up
        close all
        clear dir_subj subj freq spec h source
    end
    
    %% Get IRASA Power
    % restart timer
    tic
    
    % cycle through each subject
    for subj = 1 : n_subj

        % define subject data directory
        dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];    

        % load in raw data
        fprintf('\nloading sub-%02.0f data...\n',subj);
        load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'],'source')  

        % get irasa
        freq{1} = get_roi_irasa(source,mask_roi{1},'encoding','visual');
        freq{2} = get_roi_irasa(source,mask_roi{2},'retrieval','visual');

        % save spectral outputs
        save(sprintf('%sderivatives/sub-%02.0f/eeg/sub-%02.0f_task-corr_eeg-irasa.mat',dir_bids,subj,subj),'freq')  
        fprintf('\nsub-%02.0f complete...\ntotal time elapsed: %1.1f minutes...\nestimated time remaining: %1.1f minutes...\n',subj,round(toc/60,1),round((round(toc/60,1)/subj)*(n_subj-subj),1))

        % tidy up
        close all
        clear dir_subj subj freq spec h source
    end
end

% --- EEG DONE --------------------------------------------------------- %%

%% Correlate Similarity with Wavelet
% load in similarity index
load([dir_bids,'derivatives/group/rsa-correlation/group_task-all_fmri-si'],'rsa_vec')

% load mean bold
load([dir_bids,'derivatives/group/rsa-correlation/group_task-all_fmri-meanbold'],'mean_bold')
   
% predefine structure for correlation values
grand_freq = repmat({struct('dimord','subj_chan_freq_time',...
                    'freq',8,...
                    'label',{{'dummy'}},...
                    'time',1,...
                    'cfg',[])},[2 3]); tic
                
% cycle through each subject
for subj = 1 : n_subj
    
    % load in irasa data
    load(sprintf('%sderivatives/sub-%02.0f/eeg/sub-%02.0f_task-corr_eeg-wavelet.mat',dir_bids,subj,subj),'freq')  
         
    % cycle through each mask (plus one for forgotten ers)
    for mask = 1 : numel(mask_names)
        
        % extract all trial numbers from powspctrm
        trl_nums = freq{mask}.trialinfo(:,2+mask);

        % extract bold for these trials
        bold_tmp = squeeze(mean_bold(subj,mask,trl_nums));
        
        % adjsut trials numbers for rsa_vec
        trl_nums(trl_nums>48 & trl_nums<=96)    = trl_nums(trl_nums>48 & trl_nums<=96) - 48;
        trl_nums(trl_nums>96 & trl_nums<=144)   = trl_nums(trl_nums>96 & trl_nums<=144) - 48;
        trl_nums(trl_nums>144 & trl_nums<=192)  = trl_nums(trl_nums>144 & trl_nums<=192) - 96;
        
        % create regressor matrix
        X      = zeros(numel(trl_nums),4);
        X(:,1) = rsa_vec(subj,mask,trl_nums);
        X(:,2) = freq{mask}.powspctrm(:,:,:,1) - freq{mask}.powspctrm(:,:,:,2);
        X(:,3) = freq{mask}.powspctrm(:,:,:,2);
        X(:,4) = freq{mask}.powspctrm(:,:,:,1);
        X(:,5) = bold_tmp;
        X(:,6) = freq{mask}.trialinfo(:,2);

        % if retrieval, drop forgotten
        if mask == 2; X(freq{mask}.trialinfo(:,1)~=1,:) = []; end
       
        % create table
        tbl = array2table(X,'VariableNames',{'rsa','erd','prepow','postpow','bold','conf'});

        % --- fit ERD model --- %
        % define predictors
        preds = {'bold','conf','erd'};
        
        % run linear model
        B = fitlm(tbl,'ResponseVar','rsa',...
                      'RobustOpts','on',...
                      'PredictorVars',preds);

        % extract predictor results
        for p = 1 : numel(preds); grand_freq{mask,1}.(preds{p})(subj,1) = B.Coefficients.tStat(preds{p}); end   
        
        % --- fit pre/post model --- %
        % define predictors
        preds = {'bold','conf','postpow','prepow'};
        
        % run linear model
        B = fitlm(tbl,'ResponseVar','rsa',...
                      'RobustOpts','on',...
                      'PredictorVars',preds);

        % extract predictor results
        for p = 1 : numel(preds); grand_freq{mask,2}.(preds{p})(subj,1) = B.Coefficients.tStat(preds{p}); end   
        
        % --- fit median split model --- %
        % define predictors
        preds = {'bold','conf','postpow','prepow'};
        
        % median split variables
        tbl.prepow  = double(tbl.prepow>=median(tbl.prepow));
        tbl.postpow = double(tbl.postpow>=median(tbl.postpow));
        
        % run linear model
        B = fitlm(tbl,'ResponseVar','rsa',...
                      'RobustOpts','on',...
                      'PredictorVars',preds);

        % extract predictor results
        for p = 1 : numel(preds); grand_freq{mask,3}.(preds{p})(subj,1) = B.Coefficients.tStat(preds{p}); end     
    end
    
    % update user
    fprintf('sub-%02.0f complete...\n',subj)
end

% save data
save([dir_bids,'derivatives/group/rsa-correlation/group_task-rf_comb-freq.mat'],'grand_freq'); 

% --- Run Statistics
% set seed
rng(1) 

% define parameter and tail
tail     = cell(3,1);
tail{1}  = [0 0 -1];
tail{2}  = [0 0 -1 0];
tail{3}  = [0 0 -1 0];

% cycle through each model
for model = 1 : 3

    % predefine stat/tbl outputs
    sout = cell(numel(tail{model}),1);
    tout = cell(numel(tail{model}),1);

    % get predictors
    preds = fieldnames(grand_freq{1,model});
    preds = preds(~ismember(preds,{'dimord','freq','label','time','cfg'}));
    
    % cycle through each parameter
    for i = 1 : numel(preds)

        % predefine cell for statistics
        cfg             = [];
        cfg.tail        = [0 0] + tail{model}(i);
        cfg.parameter   = preds{i};
        [sout{i},tout{i}] = run_oneSampleT(cfg, grand_freq(:,model)); 
    end

    % combine tables
    tbl = cat(1,tout{:});

    % write cluster information
    clus_tbl = int_extractSubjs(grand_freq(:,model),preds,{'enc','ret'});
    writetable(clus_tbl,[dir_bids,'derivatives/group/eeg/group_task-rf_eegfmri-wavecluster.csv'],'Delimiter',',')
    copyfile([dir_bids,'derivatives/group/eeg/group_task-rf_eegfmri-wavecluster.csv'],[dir_repos,'data/fig3_data/group_task-rf_eegfmri-wavecluster.csv'])

    % extract data for Hotellings
    if model == 2
        hot_data(:,1:2) = [clus_tbl.enc_postpow clus_tbl.ret_postpow];
    end
    
    % add labels
    opt = repmat({'encoding','retrieval'},[1 numel(preds)])';
    var = repmat(preds',[2 1]);
    tbl = addvars(tbl,var(:),'Before','t');
    tbl = addvars(tbl,opt(:),'Before','t');
    disp(tbl)    
end

% extract data for Hotellings
save([dir_bids,'derivatives/group/eeg/group_task-rf_eegfmri-hot_data.mat'],'hot_data')

%% Correlate Similarity and IRASA
% load in similarity index
load([dir_bids,'derivatives/group/rsa-correlation/group_task-all_fmri-si'],'rsa_vec')

% load mean bold
load([dir_bids,'derivatives/group/rsa-correlation/group_task-all_fmri-meanbold'],'mean_bold')
   
% predefine structure for correlation values
grand_freq = repmat({struct('dimord','subj_chan_freq_time',...
                    'freq',8,...
                    'label',{{'dummy'}},...
                    'time',1,...
                    'cfg',[])},[2 1]); tic
                
% cycle through each subject
for subj = 1 : n_subj
    
    % load in irasa data
    load(sprintf('%sderivatives/sub-%02.0f/eeg/sub-%02.0f_task-corr_eeg-irasa.mat',dir_bids,subj,subj),'freq')  
    
    % cycle through each mask (plus one for forgotten ers)
    for mask = 1 : numel(mask_names)
        
        % extract all trial numbers from powspctrm
        trl_nums = freq{mask}.trialinfo(:,2+mask);
        
        % extract bold for these trials
        bold_tmp = squeeze(mean_bold(subj,mask,trl_nums));
        
        % fix numbers
        trl_nums(trl_nums>48 & trl_nums<=96)    = trl_nums(trl_nums>48 & trl_nums<=96) - 48;
        trl_nums(trl_nums>96 & trl_nums<=144)   = trl_nums(trl_nums>96 & trl_nums<=144) - 48;
        trl_nums(trl_nums>144 & trl_nums<=192)  = trl_nums(trl_nums>144 & trl_nums<=192) - 96;       
        
        % get alpha/beta power
        pow      = [];
        pow(:,1) = squeeze(mean(freq{mask}.powspctrm(:,1,freq{mask}.freq>=8),3)); 
        pow(:,2) = squeeze(mean(freq{mask}.powspctrm(:,2,freq{mask}.freq>=8),3)); 
        
        % predefine slope out
        frac = zeros(size(freq{mask}.fractal,1),2,2);
        
        % cycle through trials
        for trl = 1 : size(freq{mask}.fractal,1)

            % log transform fractal
            logF = log10(freq{mask}.freq');
            logP = squeeze(log10(freq{mask}.fractal(trl,:,:)))';

            % get linfit
            Afit = fitlm(logF,logP(:,1)); %post
            Bfit = fitlm(logF,logP(:,2)); %pre

            % get slope and intercept difference
            frac(trl,1,:) = Bfit.Coefficients.tStat;
            frac(trl,2,:) = Afit.Coefficients.tStat;
        end

        % create regressor matrix
        X      = zeros(numel(trl_nums),4);
        X(:,1) = rsa_vec(subj,mask,trl_nums);
        X(:,2) = pow(:,1)-pow(:,2);
        X(:,3) = pow(:,2);
        X(:,4) = frac(:,1,1)-frac(:,2,1);
        X(:,5) = frac(:,1,2)-frac(:,2,2);
        X(:,6) = frac(:,2,1);
        X(:,7) = frac(:,2,2);
        X(:,8) = bold_tmp;
        X(:,9) = freq{mask}.trialinfo(:,2);
        
        % if retrieval, drop forgotten
        if mask == 2; X(freq{mask}.trialinfo(:,1)~=1,:) = []; end
        
        % create table
        tbl = array2table(X,'VariableNames',{'rsa','postpow','prepow','postint','postslp','preint','preslp','bold','conf'});

        % define predictors
        preds = {'bold','conf','postpow','postslp','postint'};
        
        % run linear model
        B = fitlm(tbl,'ResponseVar','rsa',...
                      'PredictorVars',preds,...
                      'RobustOpts','on');
                  
        % extract predictor results
        for p = 1 : numel(preds)
                        
            % store data
            grand_freq{mask}.(preds{p})(subj,1) = B.Coefficients.tStat(preds{p});
        end 
    end
    
    % update user
    fprintf('sub-%02.0f complete...\n',subj)
end

% save data
save([dir_bids,'derivatives/group/rsa-correlation/group_task-all_comb-freq.mat'],'grand_freq'); 

% Run Statistics
% set seed
rng(1) 

% define parameter and tail
tail  = [0 0 -1 0 0 0 0 0];

% predefine stat/tbl outputs
sout = cell(numel(preds),1);
tout = cell(numel(preds),1);

% cycle through each parameter
for i = 1 : numel(preds)
    
    % predefine cell for statistics
    cfg             = [];
    cfg.tail        = [0 0] + tail(i);
    cfg.parameter   = preds{i};
    [sout{i},tout{i}] = run_oneSampleT(cfg, grand_freq(:)); 
end

% combine tables
tbl = cat(1,tout{:});

% add labels
opt = repmat({'encoding','retrieval'},[1 numel(preds)])';
var = repmat(preds,[2 1]);
tbl = addvars(tbl,var(:),'Before','t');
tbl = addvars(tbl,opt(:),'Before','t');
disp(tbl)

end

%% SUBFUNCTIONS

function tbl = int_extractSubjs(data,metric,condition)

count = 1;
output = nan(size(data{1}.(metric{1}),1),numel(metric)*numel(condition));
label = cell(numel(metric)*numel(condition),1);
for i = 1 : numel(metric)
    for j = 1 : numel(condition)
           
        % save output
        output(:,count) = data{j}.(metric{i});
        
        % define label
        label{count} = [condition{j},'_',metric{i}];
        count = count + 1;
    end
end

% create table
tbl = array2table(output,'VariableNames',label);

end

