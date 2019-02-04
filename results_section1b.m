%% RSA
% initialisation
clearvars
clc

% initalise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% define root directory
if ispc        
    dir_root = 'Y:/projects/reinstatement_fidelity/';
    dir_bids = [dir_root,'bids_data/'];
    dir_tool = 'Y:/projects/general/';
    dir_repos = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
    copylocal = true;
    
    % add subfunctions
    addpath([dir_repos,'subfunctions'])
    
else       
    dir_root = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/';
    dir_bids = [dir_root,'bids_data/'];
    dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';            
    dir_repos = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/git_clone/reinstatement_fidelity/'; % repository directory
    copylocal = false;
    
    % add subfunctions
    addpath([dir_repos,'subfunctions'])
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
scan_search = 12;                                           % searchlight radius
scan_func   = {'_3_1','_4_1','_5_1','_6_1',...
               '_8_1','_9_1','_10_1','_11_1'};              % functional scan suffix

%% Create GLM Model
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/'])
    mkdir([dir_subj,'rsa-ers/'])
    
    % delete old SPM file if it exists
    if exist([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/SPM.mat'],'file'); delete([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/SPM.mat']); end
    
    % create cell to hold events data
    events_onset  = zeros(24,8);
    count = ones(1,8);
    
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
            
            % check if encoding
            if strcmpi(tbl.operation{e},'encoding'); offset = 0;
            else; offset = 4;
            end
            
            % switch according to string
            switch tbl.stimulus{e}
                case 'WATERMILL';   idx = 1 + offset;
                case 'BIKE';        idx = 2 + offset;
                case 'UNDERWATER';  idx = 3 + offset;
                case 'FARM';        idx = 4 + offset;
                otherwise;          continue
            end
            
            % add key values
            events_onset(count(1,idx),idx) = tbl.onset(e);
            
            % sort memory
            if tbl.recalled(e) ~= 1
                events_onset(count(1,idx),idx) = NaN;
            end
            
            % get button press onset (if pressed)
            if ~isnan(tbl.rt(e)); button_onset(end+1,1) = tbl.onset(e) + tbl.rt(e); end %#ok<SAGROW>
                
            % update counter
            count(1,idx) = count(1,idx) + 1;
            
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
    else; save([dir_subj,'/rsa-ers/R.mat'],'R')   
    end
    
    % get all scans for GLM
    all_scans = get_functional_files([dir_root,'bids_data/derivatives/',subj_handle,'/'],'ua',copylocal);
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
        matlabbatch{1}.spm.stats.fmri_spec.dir              = {[dir_subj,'/rsa-ers']};
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {[dir_subj,'/rsa-ers/R.mat']}; 
    end
    
    % cycle through and define each condition
    for trl = 1 : size(events_onset,2)
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(trl).name        = ['trl',sprintf('%03.0f',trl)];
        matlabbatch{1}.spm.stats.fmri_spec.sess.cond(trl).onset       = events_onset(~isnan(events_onset(:,trl)),trl);
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
    else; matlabbatch{2}.spm.stats.fmri_est.spmmat = {[dir_subj,'/rsa-ers/SPM.mat']};
    end
    
    % run batch
    spm_jobman('run',matlabbatch)
    clear matlabbatch all_scans bad_scans subj_handle button_onset events_onset trl    
        
    % copy and delete local files if they exist
    if copylocal
        copyfile('C:/tmp_mri/spm/R.mat',[dir_subj,'/rsa-ers/R.mat'])
        copyfile('C:/tmp_mri/spm/SPM.mat',[dir_subj,'/rsa-ers/SPM.mat'])
        copyfile('C:/tmp_mri/spm/mask.nii',[dir_subj,'/rsa-ers/mask.nii'])
        rmdir('C:/tmp_mri/','s'); 
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
    
    % predefine matrix for mask data
    maskImg = zeros(1, prod(scan_fov));
    
    % load mask and add to matrix    
    nii = load_untouch_nii([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/mask.nii']);
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
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/'])
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/',subj_handle,'_task-rf_rsa-fullVolume.mat'],'scanVec')
    
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
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/',subj_handle,'_task-rf_rsa-maskedVolume.mat'],'patterns')
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/',subj_handle,'_task-rf_rsa-mask.mat'],'mask_idx')
    
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
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/',subj_handle,'_task-rf_rsa-maskedVolume.mat'])
    
    % get mean pattern across all trials and replicate matrix
    meanPattern = repmat(mean(patterns,1),[size(patterns,1), 1]);

    % subtract mean pattern from data
    patterns = patterns - meanPattern;
    
    % save patterns
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/',subj_handle,'_task-rf_rsa-maskedDemeanedVolume.mat'],'patterns')
    
    % clean up
    clear subjHandle meanPattern patterns
end

%% Prepare GLM and Data for Searchlight Analysis
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    
    % load pattern data
    load([dir_subj,'/rsa-ers/',subj_handle,'_task-rf_rsa-maskedDemeanedVolume.mat'])
    
    % load SPM.mat
    load([dir_subj,'/rsa-ers/SPM.mat'])
    
    % create table to record stimulus detail
    stim_details  = array2table(zeros(n_trials*2,2),'VariableNames', {'encoding','modality'});
    stim_count = 1;
    
    % create table to record scan details
    scan_details  = array2table(zeros(n_volumes*8,2),'VariableNames', {'encoding','modality'});
    scan_count = 1;  
    
    % cycle through each run
    for run = 1 : n_runs

        % load event table
        tbl = readtable([dir_root,'bids_data/',subj_handle,'/func/',subj_handle,'_task-rf_run-',num2str(run),'_events.tsv'],'FileType','text','Delimiter','\t');
        
        % check block types
        block_encoding  = double(any(ismember(tbl.operation,'encoding')));
        block_visual    = double(any(ismember(tbl.modality,'Visual')));
        
        % cycle through every event
        for e = 1 : size(tbl,1)
            
            % if an event
            if strcmpi(tbl.trial_type(e),'Stimulus Onset')
            
                % add key values
                stim_details.encoding(stim_count) = strcmpi(tbl.operation(e),'encoding');
                stim_details.modality(stim_count) = strcmpi(tbl.modality(e),'Visual');
                
                % get stimulus values
                switch tbl.stimulus{e}
                    case 'WATERMILL';   stim_details.stimulus(stim_count,1) = 1;
                    case 'UNDERWATER';  stim_details.stimulus(stim_count,1) = 2;
                    case 'BIKE';        stim_details.stimulus(stim_count,1) = 3;
                    case 'FARM';        stim_details.stimulus(stim_count,1) = 4;
                end
                
                % change value if retrieval
                if stim_details.encoding(stim_count) ~= 1 && stim_details.modality(stim_count) == 1
                    stim_details.stimulus(stim_count,1) = stim_details.stimulus(stim_count,1) + 4;
                end
                
                % update counter
                stim_count = stim_count + 1; 
                
            % else if a volume
            elseif strncmpi(tbl.trial_type(e),'Volume',6)
                
                % add key values
                scan_details.encoding(scan_count) = block_encoding;
                scan_details.modality(scan_count) = block_visual;

                % update counter
                scan_count = scan_count + 1;                 
            end          
        end
    end
    
    % kick out first three/last five scans of each run
    scan_details([1:3 251:258 506:513 761:768 1016:1023 1271:1278 1526:1533 1781:1788 2036:2040],:) = [];
    
    % clean up
    clear run tbl e stim_count
    
    % get design matrix (X) and split into two groups (Xa and Xb)
    X.raw = SPM.xX.X;
        
    % remove scans/regressors that are not visual (to
    % computationally demanding to anything more than this)
    X.raw = X.raw(scan_details.modality==1,:);
    
    % split GLM into two groups (train [encoding] and test [retrieval] data)
    X.a = X.raw(1:size(X.raw,1)/2,:);
    X.b = X.raw((size(X.raw,1)/2)+1:end,:);
    
    % --- prepare patterns --- %
    % get activation matrix (Y) and 
    Y.raw = patterns(scan_details.modality==1,:);
    
    % split into two groups (Ya and Yb)
    Y.a = Y.raw(1:size(Y.raw,1)/2,:);
    Y.b = Y.raw(size(Y.raw,1)/2+1:end,:);
        
    % --- remove singular dimensions --- %
    % find zero-value rows
    X.a_badRow = all(X.a(:,1:8)==0,2);
    X.b_badRow = all(X.b(:,1:8)==0,2);
    
    % remove zero-value rows
    X.a(X.a_badRow,:) = [];
    X.b(X.a_badRow,:) = [];
    
    % remove zero-value rows
    Y.a(X.a_badRow,:) = [];
    Y.b(X.b_badRow,:) = [];
        
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
    
    % save data
    save([dir_subj,'/rsa-ers/',subj_handle,'_task-rf_rsa-FormattedVolume.mat'],'X','Y')
    
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

% clean up
clear x y z voxRadInSearchLight

%% Run Searchlight LDt Analysis
tic

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % load data
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/',subj_handle,'_task-rf_rsa-FormattedVolume.mat'])
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/',subj_handle,'_task-rf_rsa-mask.mat'])
    
    % predefine model rdm
    model_rdm{1}        	= nan(8,8);

    % set items that belong to the different category to 1
    model_rdm{1}(5:8,1:4)	= 1; 

    % set items that belong to a same category to -1
    model_rdm{1}(5,1)       = -1;
    model_rdm{1}(6,2)       = -1;
    model_rdm{1}(7,3)       = -1;
    model_rdm{1}(8,4)       = -1;
    
    % set diagonal to NaN
    model_rdm{1}(logical(eye(size(model_rdm{1},1)))) = 0; 

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
    RDM_ldt = rsa.stat.fisherDiscrTRDM_searchlight(X.a,Y.a,X.b,Y.b,1:8,sl_vox,model_rdm,false);

    % save raw output
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/',subj_handle,'_task-rf_rsa-rawRDM.mat'],'RDM_ldt')
    
    % clean up
    clear X Y sl_vox
    
    % add z-value to rdmBrain
    rdmBrain(M(goodSL)) = mean(cat(2,RDM_ldt.ats,RDM_ldt.bts),2);

    % load template nifti
    filename = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/mask.nii'];        
    V = load_untouch_nii(filename);

    % change filename, datatype, and image
    V.fileprefix = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/',subj_handle,'_task-rf_rsa-searchlight'];
    V.hdr.dime.datatype = 64;
    V.img = rdmBrain;

    % save z-value brain
    save_untouch_nii(V,[V.fileprefix,'.nii']);  
    
    % update command line
    tElapse  = toc;
    tPerLoop = tElapse / subj;
    loopsRem = n_subj - subj;
    timeRem  = (tPerLoop * loopsRem)/3600;
    fprintf('\nSubject %1.0f of %1.0f complete...\nApproximate time remaining: %1.0f hours...\n',subj,n_subj,timeRem)
           
    % tidy
    clear V filename rdmBrain avgZ mN i goodSL M mask_idx RDM_ldt subjHandle
end

% clean up
clear subj tElapse tPerLoop loopsRem timeRem

%% Normalise and Smooth
% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % normalise functional
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb      = [-78 -112 -70; 78 76 85];    
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox     = [3 3 4];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp  = 4;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def         = {[dir_root,'bids_data/derivatives/',subj_handle,'/anat/y_',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample    = {[dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/',subj_handle,'_task-rf_rsa-searchlight.nii,1']};

    % smooth
    matlabbatch{2}.spm.spatial.smooth.fwhm                      = [8 8 8];
    matlabbatch{2}.spm.spatial.smooth.dtype                     = 0;
    matlabbatch{2}.spm.spatial.smooth.im                        = 1;
    matlabbatch{2}.spm.spatial.smooth.prefix                    = 's';
    matlabbatch{2}.spm.spatial.smooth.data                      = {[dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/w',subj_handle,'_task-rf_rsa-searchlight.nii,1']};
    
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
    rMapFiles{subj,1}  = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/sw',subj_handle,'_task-rf_rsa-searchlight.nii,1'];
end

% create second-level glms
matlabbatch{1}.spm.stats.factorial_design.dir                       = {[dir_root,'bids_data/derivatives/group/rsa-ers']};
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
matlabbatch{1}.spm.stats.fmri_est.spmmat            = {[dir_root,'bids_data/derivatives/group/rsa-ers/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals   = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1;

spm_jobman('run',matlabbatch)
clear matlabbatch

% define contrasts
matlabbatch{1}.spm.stats.con.spmmat(1)                  = {[dir_root,'bids_data/derivatives/group/rsa-ers/SPM.mat']};   
matlabbatch{1}.spm.stats.con.delete                     = 0;    
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name       = 'within>between';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec     = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none';  

spm_jobman('run',matlabbatch)
clear matlabbatch

%% Extract Metrics for Visualisation
% combine visual cluster and save
combine_spm_cluster([dir_root,'bids_data/derivatives/group/rsa-ers/'])

% load SPM details
load([dir_root,'bids_data/derivatives/group/rsa-ers/SPM.mat'])

% extract subject values for each cluster
[betas,d] = extract_sample_points([dir_root,'bids_data/derivatives/group/rsa-ers/'],SPM);

% save betas as table
tbl = array2table(betas','VariableNames',{'LeftHemi','RightHemi'});
writetable(tbl,[dir_repos,'data/ers_betas.csv'],'Delimiter',',')

% save effect size as table
tbl = array2table(d','VariableNames',{'LeftHemi','RightHemi'});
writetable(tbl,[dir_repos,'data/ers_cohensD.csv'],'Delimiter',',')
