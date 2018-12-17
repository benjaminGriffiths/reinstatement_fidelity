%% RSA
% initialisation
clearvars
clc

% initalise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% define root directory
if ispc;        
    dir_root = 'Y:/projects/reinstatement_fidelity/';
    dir_bids = [dir_root,'bids_data/'];
    dir_tool = 'Y:/projects/general/';
    dir_repos = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
    
    % add subfunctions
    addpath([dir_repos,'subfunctions'])
    
else       
    dir_root = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/';
    dir_bids = [dir_root,'bids_data/'];
    dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';            
    dir_repos = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/scripts/'; % repository directory
end

% add RSA toolbox to path
addpath(genpath([dir_tool,'rsatoolbox-develop']))

%% Define Key Parameters
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

%% Prepare Masks
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/',subj_handle,'/'];
    
    % make directory for masks
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/masks/'])
    
    % prepare deformation batch
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def        = {[dir_root,'bids_data/derivatives/',subj_handle,'/anat/iy_',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.util.defs.comp{1}.inv.space              = {[dir_root,'bids_data/',subj_handle,'/anat/',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames             = {[dir_root,'bids_data/sourcedata/masks/whole_brain.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.weight             = {''};
    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr    = {[dir_root,'bids_data/derivatives/',subj_handle,'/masks/']};
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.file           = {[dir_root,'bids_data/derivatives/',subj_handle,'/func/meanua',subj_handle,'_task-rf_run-1_bold.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.preserve           = 0;
    matlabbatch{1}.spm.util.defs.out{1}.push.fwhm               = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.push.prefix             = '';
    
    % run
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
    % move file to subject directory
    movefile([dir_root,'bids_data/derivatives/',subj_handle,'/masks/wwhole_brain.nii'],...
        [dir_root,'bids_data/derivatives/',subj_handle,'/masks/whole_brain.nii'])
    
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
    nii = load_untouch_nii([dir_root,'bids_data/derivatives/',subj_handle,'/masks/whole_brain.nii']);
    maskImg(1,:) = reshape(nii.img,1,[]);

    % predefine matrix for functional data
    scanVec = zeros((n_volumes-8).*numel(scan_func),numel(nii.img));
    
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
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/'])
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-fullVolume.mat'],'scanVec')
    
    % apply mask to scans
    fprintf('Masking data...\n')
    
    % define patterns for mask
    patterns = scanVec;

    % mask functional data (set all zero-element voxels in mask to zero)
    patterns(:,maskImg(1,:)==0) = 0;

    % locate "dead voxels" (functional data within mask that has zero values; e.g. where mask captures skull)
    deadIdx = any(patterns == 0);

    % set these voxels to zero across all scans
    patterns(:,deadIdx) = 0;
    
    % get a boolean vector of non-zero voxels in patterns matrix
    mask_idx = all(patterns~=0);
    
    % remove zero elements from patterns
    patterns(:,any(patterns == 0)) = [];
    
    % save
    fprintf('Saving masked volumes...\n')
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-maskedVolume.mat'],'patterns')
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-mask.mat'],'mask_idx')
    
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
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-maskedVolume.mat'])
    
    % get mean pattern across all trials and replicate matrix
    meanPattern = repmat(mean(patterns,1),[size(patterns,1), 1]);

    % subtract mean pattern from data
    patterns = patterns - meanPattern;
    
    % save patterns
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-maskedDemeanedVolume.mat'],'patterns')
    
    % clean up
    clear subjHandle meanPattern patterns
end

%% Create GLM Model
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/',subj_handle,'/'];
    
    % delete old SPM file if it exists
    if exist([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/SPM.mat'],'file'); delete([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/SPM.mat']); end
    
    % create cell to hold events data
    events_onset  = zeros(n_trials*2,1);
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
            
            % add key values
            events_onset(count) = tbl.onset(e);
            
            % get button press onset (if pressed)
            if ~isnan(tbl.rt(e)); button_onset(end+1,1) = tbl.onset(e) + tbl.rt(e); end %#ok<SAGROW>
                
            % update counter
            count = count + 1;
            
        end
    end
    
    % tidy up
    clear count e tbl LoS run
       
    % convert event onsets and button presses from EEG samples to seconds
    events_onset = events_onset ./ EEG_sample;
    button_onset = button_onset ./ EEG_sample;    
    
    % load movement parameters
    R = load([dir_root,'bids_data/derivatives/',subj_handle,'/func/rp_a',subj_handle,'_task-rf_run-1_bold.txt']);
    
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
    R(1+(n_volumes_adj*0):n_volumes_adj*1,15) = zeros(1,n_volumes_adj);
    R(1+(n_volumes_adj*2):n_volumes_adj*3,16) = zeros(1,n_volumes_adj);
    R(1+(n_volumes_adj*4):n_volumes_adj*5,17) = zeros(1,n_volumes_adj);
    R(1+(n_volumes_adj*6):n_volumes_adj*7,18) = zeros(1,n_volumes_adj);
    R(1+(n_volumes_adj*1):n_volumes_adj*2,19) = zeros(1,n_volumes_adj);
    R(1+(n_volumes_adj*3):n_volumes_adj*4,20) = zeros(1,n_volumes_adj);
    R(1+(n_volumes_adj*5):n_volumes_adj*6,21) = zeros(1,n_volumes_adj);
    R(1+(n_volumes_adj*7):n_volumes_adj*8,21) = zeros(1,n_volumes_adj); 
    
    % save nuisance regressors
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/R.mat'],'R')            
    clear R n_volumes_adj run
    
    % get all scans for GLM
    all_scans = get_functional_files([dir_root,'bids_data/derivatives/',subj_handle,'/'],'ua');
    all_scans = all_scans{1};
    
    % remove bad scans
    all_scans(bad_scans) = [];
    
    % define parameters for GLM
    matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
    matlabbatch{1}.spm.stats.fmri_spec.dir              = {[dir_root,'bids_data/derivatives/',subj_handle,'/rsa']};
    matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi       = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress     = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf         = 128;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans       = all_scans;
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {[dir_root,'bids_data/derivatives/',subj_handle,'/rsa/R.mat']};   
    matlabbatch{1}.spm.stats.fmri_spec.timing.units     = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT        = 2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t    = 32;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0   = 16;  
                   
    % cycle through and define each condition
    for trl = 1 : numel(events_onset)
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
    
    % specify model
    spm_jobman('run',matlabbatch)
    clear matlabbatch all_scans bad_scans subj_handle button_onset events_onset trl
    
end

%% Prepare GLM and Data for Searchlight Analysis
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    
    % load pattern data
    load([dir_subj,'/rsa/',subj_handle,'_task-rf_rsa-maskedDemeanedVolume.mat'])
    
    % load SPM.mat
    load([dir_subj,'/rsa/SPM.mat'])
    
    % create table to record stimulus detail
    stim_details  = array2table(zeros(n_trials*2,2),'VariableNames', {'encoding','remembered'});
    stim_details.word = cell(n_trials*2,1);
    stim_details.dynamic = cell(n_trials*2,1);
    count = 1;
    
    % cycle through each run
    for run = 1 : n_runs

        % load event table
        tbl = readtable([dir_root,'bids_data/',subj_handle,'/func/',subj_handle,'_task-rf_run-',num2str(run),'_events.tsv'],'FileType','text','Delimiter','\t');
        
        % cut table to stimulus onset
        tbl = tbl(ismember(tbl.trial_type,'Stimulus Onset'),:);
         
        % cycle through every event
        for e = 1 : size(tbl,1)
            
            % add key values
            stim_details.word(count)       = tbl.word(e);
            stim_details.dynamic(count)    = tbl.stimulus(e);
            stim_details.encoding(count)   = strcmpi(tbl.operation(e),'encoding');
            stim_details.remembered(count) = tbl.recalled(e);
                       
            % update counter
            count = count + 1;            
        end
    end
    
    % clean up
    clear run tbl e count
    
    % get design matrix (X) and split into two groups (Xa and Xb)
    X.raw = SPM.xX.X;
    
    % get matrix of trials by stimulus content at encoding
    trl_by_con{1} = find(strcmpi(stim_details.dynamic,'bike') & stim_details.encoding == 1);
    trl_by_con{2} = find(strcmpi(stim_details.dynamic,'farm') & stim_details.encoding == 1);
    trl_by_con{3} = find(strcmpi(stim_details.dynamic,'underwater') & stim_details.encoding == 1);
    trl_by_con{4} = find(strcmpi(stim_details.dynamic,'watermill') & stim_details.encoding == 1);
        
    % get matrix of trials by stimulus content at retrieval
    trl_by_con{5} = find(strcmpi(stim_details.dynamic,'bike') & stim_details.encoding == 0);
    trl_by_con{6} = find(strcmpi(stim_details.dynamic,'farm') & stim_details.encoding == 0);
    trl_by_con{7} = find(strcmpi(stim_details.dynamic,'underwater') & stim_details.encoding == 0);
    trl_by_con{8} = find(strcmpi(stim_details.dynamic,'watermill') & stim_details.encoding == 0);
    
    % get GLM for nuisance regressors
    X.n = X.raw(:,385:end);
    
    % get condition average GLM
    for i = 1 : numel(trl_by_con)
        X.c(:,i) = sum(X.raw(:,trl_by_con{i}),2);
    end
    
    % clean up
    clear i trl_by_con SPM stim_details
    
    % split GLM into two groups (train and test data) [excluding nuisance regressors]
    X.a = X.c(1:size(X.raw)/2,:);
    X.b = X.c((size(X.raw)/2)+1:end,:);
    
    % split nuisance regressors into two groups
    X.n_a = X.n(1:size(X.raw)/2,:);
    X.n_b = X.n((size(X.raw)/2)+1:end,:);
    
    % find zero-value rows
    X.a_badRow = all(X.a==0,2);
    X.b_badRow = all(X.b==0,2);
    
    % remove zero-value rows
    X.a(X.a_badRow,:)   = [];
    X.b(X.b_badRow,:)   = [];    
    X.n_a(X.a_badRow,:) = [];
    X.n_b(X.b_badRow,:) = [];
    
    % add nuisance regressors
    X.a(:,end+1:end+size(X.n_a,2)) = X.n_a;
    X.b(:,end+1:end+size(X.n_b,2)) = X.n_b;
        
    % cycle through A/B avg/trl combinations
    fn = {'a','b'};
    for i = 1 : 2
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
    X.a(:,X.a_badCol) = [];
    X.b(:,X.b_badCol) = [];
    
    % get activation matrix (Y) and split into two groups (Ya and Yb)
    Y.raw = patterns;
    Y.a = Y.raw(1:size(Y.raw,1)/2,:);
    Y.b = Y.raw((size(Y.raw,1)/2)+1:end,:);

    % remove zero-value rows
    Y.a(X.a_badRow,:) = [];
    Y.b(X.b_badRow,:) = [];

    % save data
    save([dir_subj,'/rsa/',subj_handle,'_task-rf_rsa-FormattedVolume.mat'],'X','Y')
    
    clear X Y i j fn patterns subj_handle dir_subj
    
    fprintf('Subject %02.0f of %02.0f prepared...\n',subj,n_subj)
end

%% Define Searchlight Characteristics
% get the radius of searchlight in voxels
voxRadInSearchLight = scan_search ./ scan_vox;

% define distance from searchlight centre to perimeter in voxels
dist2Perimeter = ceil(voxRadInSearchLight);

% create boolean searchlight sphere
[x,y,z] = meshgrid(-dist2Perimeter(1):dist2Perimeter(1),-dist2Perimeter(2):dist2Perimeter(2),-dist2Perimeter(3):dist2Perimeter(3));
sphere  = ((x*scan_vox(1)).^2+(y*scan_vox(2)).^2+(z*scan_vox(3)).^2)<=(scan_search^2);

% predefine model rdm
model_rdm{1} = nan(8,8);

% set items that belong to the different category to 1
model_rdm{1}(5:8,1:4) = 1;

% set items that belong to a same category to -1
model_rdm{1}(5,1) = -1;
model_rdm{1}(6,2) = -1;
model_rdm{1}(7,3) = -1;
model_rdm{1}(8,4) = -1;

% set diagonal to NaN
model_rdm{1}(logical(eye(size(model_rdm{1},1)))) = 0; 

% clean up
clear x y z voxRadInSearchLight

%% Run Searchlight LDt Analysis
tic

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % load data
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-FormattedVolume.mat'])
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-mask.mat'])
    
    % define empty brain map
    rdmBrain = zeros(scan_fov);
    
    % predefine voxels to be analysed
    sl_vox = cell(sum(mask_idx),1);
    
    % predefine vector of searchlights to be used
    goodSL = true(size(sl_vox));
    
    % get index of all voxels within mask
    M = find(mask_idx);
    
    % update command line
    fprintf('Defining searchlights...\n')
    
    % run LDt analysis for every voxel
    for vox = 1 : numel(M)
        
        % get subscript co-ordinates of searchlight centre
        [x,y,z] = ind2sub(scan_fov,M(vox));
        
        % define box which houses spherical searchlight
        tmpX = x-dist2Perimeter(1):x+dist2Perimeter(1);
        tmpY = y-dist2Perimeter(2):y+dist2Perimeter(2);
        tmpZ = z-dist2Perimeter(3):z+dist2Perimeter(3);
        
        % remove co-ordinates less than or equal to zero, and greater than FOV
        tmpX(tmpX<=0) = []; tmpX(tmpX>scan_fov(1)) = [];
        tmpY(tmpY<=0) = []; tmpY(tmpY>scan_fov(2)) = [];
        tmpZ(tmpZ<=0) = []; tmpZ(tmpZ>scan_fov(3)) = [];
        
        % get indices of every voxel in sphere
        voxIdx = [];
        for xi = 1 : numel(tmpX)
            for yi = 1 : numel(tmpY)
                for zi = 1 : numel(tmpZ)
                    if sphere(xi,yi,zi)
                        voxIdx(end+1) = sub2ind(scan_fov,tmpX(xi),tmpY(yi),tmpZ(zi)); %#ok<SAGROW>
                    end
                end
            end
        end
        
        % remove those indices which are not included in mask
        voxIdx(~ismember(voxIdx,M)) = []; %#ok<SAGROW>
               
        % convert whole brain indices to mask based indices
        sl_vox{vox} = find(ismember(M,voxIdx));
        
        % caluculate number of voxels in searchlight
        percentIncluded = numel(sl_vox{vox})/sum(reshape(sphere==1,1,[]));
        if percentIncluded < 0.6 % if less than 60% of max
            goodSL(vox) = false;
        end        
    end       
    
    % clean up
    clear percentIncluded voxIdx xi yi zi tmpX tmpY tmpZ vox

    % remove searchlights with less than threshold number of voxels
    sl_vox(goodSL==false) = [];
    
    % run LDt analysis
    RDM_ldt = rsa.stat.fisherDiscrTRDM_searchlight(X.a,Y.a,X.b,Y.b,1:8,sl_vox,model_rdm);

    % clean up
    clear X Y sl_vox
    
    % get mean corrcoef
    avgZ = mean(cat(2,RDM_ldt.ats(:,1),RDM_ldt.bts(:,1)),2);

    % add z-value to rdmBrain
    rdmBrain(M(goodSL)) = avgZ;

    % load template nifti
    filename = [dir_root,'bids_data/derivatives/',subj_handle,'/func/meanua',subj_handle,'_task-rf_run-1_bold.nii'];        
    V = load_untouch_nii(filename);

    % change filename, datatype, and image
    V.fileprefix = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-searchlight'];
    V.hdr.dime.datatype = 64;
    V.img = rdmBrain;

    % save z-value brain
    save_untouch_nii(V,[V.fileprefix,'.nii']);  
    
    % update command line
    tElapse  = toc;
    tPerLoop = tElapse / subj;
    loopsRem = numel(n_subj) - subj;
    timeRem  = (tPerLoop * loopsRem)*3600;
    fprintf('\nSubject %1.0f of %1.0f complete...\nApproximate time remaining: %1.0f hours...\n',subj,n_subj,timeRem)
           
    % tidy
    clear V filename rdmBrain avgZ mN i goodSL M mask_idx RDM_ldt subjHandle
end

% clean up
clear subj tElapse tPerLoop loopsRem timeRem

%% Normalise and Smooth
% define images to analyse
mN = {'Visual','Auditory'};

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % normalise functional
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb      = [-78 -112 -70; 78 76 85];    
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox     = [3 3 4];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp  = 4;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def         = {[dir_root,'bids_data/derivatives/',subj_handle,'/anat/y_',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample    = {[dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-searchlight',mN{1},'.nii,1'];
                                                                   [dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-searchlight',mN{2},'.nii,1']};

    % smooth
    matlabbatch{2}.spm.spatial.smooth.fwhm                      = [8 8 8];
    matlabbatch{2}.spm.spatial.smooth.dtype                     = 0;
    matlabbatch{2}.spm.spatial.smooth.im                        = 0;
    matlabbatch{2}.spm.spatial.smooth.prefix                    = 's';
    matlabbatch{2}.spm.spatial.smooth.data                      = {[dir_root,'bids_data/derivatives/',subj_handle,'/rsa/w',subj_handle,'_task-rf_rsa-searchlight',mN{1},'.nii,1'];
                                                                   [dir_root,'bids_data/derivatives/',subj_handle,'/rsa/w',subj_handle,'_task-rf_rsa-searchlight',mN{2},'.nii,1']};
    
    % run batch
    spm_jobman('run',matlabbatch)
    clear matlabbatch
end

%% Run Second-Level Statistics
% predefine cell for searchlight image files
rMapFiles = cell(n_subj,2);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % get searchlight images
    rMapFiles{subj,1}  = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa/sw',subj_handle,'_task-rf_rsa-searchlightVisual.nii,1'];
    rMapFiles{subj,2}  = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa/sw',subj_handle,'_task-rf_rsa-searchlightAuditory.nii,1'];
end

% create second-level glm
matlabbatch{1}.spm.stats.factorial_design.dir                       = {[dir_root,'bids_data/derivatives/group/rsa/visual/']};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans              = rMapFiles(:,1);
matlabbatch{1}.spm.stats.factorial_design.cov                       = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov                 = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none        = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im                = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em                = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm           = 1;
spm_jobman('run',matlabbatch)
clear matlabbatch subjHandle subj

% estimate model
matlabbatch{1}.spm.stats.fmri_est.spmmat            = {[dir_root,'bids_data/derivatives/group/rsa/visual/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals   = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1;

spm_jobman('run',matlabbatch)
clear matlabbatch

% define contrasts
matlabbatch{1}.spm.stats.con.spmmat(1)                  = {[dir_root,'bids_data/derivatives/group/rsa/visual/SPM.mat']};   
matlabbatch{1}.spm.stats.con.delete                     = 0;    
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name       = 'within>between';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec     = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none';  

spm_jobman('run',matlabbatch)
clear matlabbatch

% create second-level glm
matlabbatch{1}.spm.stats.factorial_design.dir                       = {[dir_root,'bids_data/derivatives/group/rsa/auditory/']};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans              = rMapFiles(:,2);
matlabbatch{1}.spm.stats.factorial_design.cov                       = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov                 = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none        = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im                = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em                = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm           = 1;
spm_jobman('run',matlabbatch)
clear matlabbatch

% estimate model
matlabbatch{1}.spm.stats.fmri_est.spmmat            = {[dir_root,'bids_data/derivatives/group/rsa/auditory/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals   = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1;

spm_jobman('run',matlabbatch)
clear matlabbatch

% define contrasts
matlabbatch{1}.spm.stats.con.spmmat(1)                  = {[dir_root,'bids_data/derivatives/group/rsa/auditory/SPM.mat']};   
matlabbatch{1}.spm.stats.con.delete                     = 0;    
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name       = 'within>between';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec     = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none';  

spm_jobman('run',matlabbatch)
clear matlabbatch rMapFiles

%% Extract Data from Cluster
% load SPM details
load([dir_root,'bids_data/derivatives/group/rsa/visual/SPM.mat'])

% cycle through each cluster
for c = 1 : 2

    % load in cluster
    nii = load_nii([dir_root,'bids_data/derivatives/group/rsa/visual/cluster',num2str(c),'.nii']);

    % get 3 x m representation of cluster
    roi = [];
    [roi(1,:),roi(2,:),roi(3,:)] = ind2sub(size(nii.img),find(nii.img>0));

    % get average of cluster
    betas(c,:) = mean(spm_get_data(SPM.xY.P,roi),2); 
end

% rename each beta and add to table
leftFusiform    = betas(1,:)';
rightFusiform   = betas(2,:)';

% add to table
data = table(leftFusiform,rightFusiform);

% calculate cohens d
for b = 1 : 2
    
    % calculate cohen's dz
    X = betas(b,:);
    numer = mean(X);
    denom = sqrt(sum((X - mean(X)).^2) ./ (numel(X)-1));
    d(b,1) = numer ./ denom;
end

% export
writetable(data,[dir_repos,'data/searchlight_betas.csv'],'Delimiter',',')

