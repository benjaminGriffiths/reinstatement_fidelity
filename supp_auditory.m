function supp_auditory

% initalise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% define root directory
if ispc     
    dir_root = 'Y:/projects/reinstatement_fidelity/';
    dir_tool = 'Y:/projects/general/';
    dir_repos = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
    
    % add subfunctions
    addpath([dir_repos,'subfunctions'])
    
else       
    dir_root = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/';
    dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';            
    dir_repos = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/scripts/'; % repository directory
    
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
scan_search = 10;                                           % searchlight radius

%% Create GLM Model
% start timer
tic

% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    mkdir([dir_subj,'rsa-percept-audio/'])
    
    % delete old SPM file if it exists
    if exist([dir_subj,'rsa-percept-audio/SPM.mat'],'file'); delete([dir_subj,'rsa-percept-audio/SPM.mat']); end
    
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
            if strcmpi(tbl.operation{e},'encoding') && strcmpi(tbl.modality{e},'auditory')

                % add key values
                events_onset(count,1) = tbl.onset(e);
                count = count + 1;
            end

            % get button press onset (if pressed)
            if ~isnan(tbl.rt(e)); button_onset(end+1,1) = tbl.onset(e) + tbl.rt(e); end       %#ok<AGROW>
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
    save([dir_subj,'/rsa-percept-audio/R.mat'],'R')            
    clear R n_volumes_adj run
    
    % get all scans for GLM
    all_scans = get_functional_files(dir_subj,'ua');
    all_scans = all_scans{1};
    
    % remove bad scans
    all_scans(bad_scans) = [];
    
    % define parameters for GLM
    matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
    matlabbatch{1}.spm.stats.fmri_spec.dir              = {[dir_subj,'/rsa-percept-audio']};
    matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi       = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress     = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf         = 128;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans       = all_scans;
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {[dir_subj,'/rsa-percept-audio/R.mat']};   
    matlabbatch{1}.spm.stats.fmri_spec.timing.units     = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT        = 2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t    = 32;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0   = 16;  
                   
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
    matlabbatch{2}.spm.stats.fmri_est.spmmat                            = {[dir_subj,'/rsa-percept-audio/SPM.mat']};

    % run batch
    spm_jobman('run',matlabbatch)
    clear matlabbatch all_scans bad_scans subj_handle button_onset events_onset trl    
    
    % update user
    fprintf('\nsub-%02.0f complete, time elapsed: %1.1f hours...\n\n',subj,toc/3600)
end

%% Read Data
% cycle through each subject
for subj = 1 : n_subj
    
    % update command line
    fprintf('\n--- Working on Subject %d ---------\n',subj)
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % predefine matrix for mask data
    maskImg = zeros(1, prod(scan_fov));
    
    % load mask and add to matrix    
    nii = load_untouch_nii([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-percept-audio/mask.nii']);
    maskImg(1,:) = reshape(nii.img,1,[]);

    % predefine matrix for functional data
    scanVec = zeros((n_volumes-8).*n_runs,numel(nii.img));
    
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
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-percept-audio/',subj_handle,'_task-rf_rsa-fullVolume.mat'],'scanVec')
    
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
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-percept-audio/',subj_handle,'_task-rf_rsa-maskedVolume.mat'],'patterns')
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-percept-audio/',subj_handle,'_task-rf_rsa-mask.mat'],'mask_idx')
    
    % clear excess variables
    clear scanVec nii deadIdx patterns mask_idx subjHandle maskImg
end

%% Mean Pattern Subtraction
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % load pattern data
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-percept-audio/',subj_handle,'_task-rf_rsa-maskedVolume.mat'],'patterns')
    
    % get mean pattern across all trials and replicate matrix
    meanPattern = repmat(mean(patterns,1),[size(patterns,1), 1]);

    % subtract mean pattern from data
    patterns = patterns - meanPattern;
    
    % save patterns
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-percept-audio/',subj_handle,'_task-rf_rsa-maskedDemeanedVolume.mat'],'patterns')
    
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
    load([dir_subj,'/rsa-percept-audio/',subj_handle,'_task-rf_rsa-maskedDemeanedVolume.mat'],'patterns')
    
    % load SPM.mat
    load([dir_subj,'/rsa-percept-audio/SPM.mat'],'SPM')
    
    % create table to record stimulus detail
    scan_details  = array2table(zeros(n_volumes*8,2),'VariableNames', {'encoding','modality'});
    stim_details  = array2table(zeros(n_trials*2,3),'VariableNames', {'encoding','modality','stimulus'});
    scan_count = 1;    
    stim_count = 1;
    
    % cycle through each run
    for run = 1 : n_runs

        % load event table
        tbl = readtable([dir_root,'bids_data/',subj_handle,'/func/',subj_handle,'_task-rf_run-',num2str(run),'_events.tsv'],'FileType','text','Delimiter','\t');
        
        % check block types
        block_encoding  = double(any(ismember(tbl.operation,'encoding')));
        block_visual    = double(any(ismember(tbl.modality,'Auditory')));
        
        % cycle through every event
        for e = 1 : size(tbl,1)
            
            % if an event
            if strcmpi(tbl.trial_type(e),'Stimulus Onset')
            
                % add key values
                stim_details.encoding(stim_count) = strcmpi(tbl.operation(e),'encoding');
                stim_details.modality(stim_count) = strcmpi(tbl.modality(e),'auditory');
                
                % get stimulus values
                switch tbl.stimulus{e}
                    case 'TRUMPET';     stim_details.stimulus(stim_count,1) = 1;
                    case 'PIANO';       stim_details.stimulus(stim_count,1) = 2;
                    case 'GUITAR';      stim_details.stimulus(stim_count,1) = 3;
                    case 'ACCORDIAN';   stim_details.stimulus(stim_count,1) = 4;
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
    clear run tbl e scan_count stim_count
    
    % --- prepare GLM --- %
    % get design matrix (X) and split into two groups (Xa and Xb)
    X.raw = SPM.xX.X;
    
    % remove scans/regressors that are not visual (to
    % computationally demanding to anything more than this)
    X.raw = X.raw(scan_details.modality==1&scan_details.encoding==1,:);
    
    % split GLM into two groups (train [encoding] and test [retrieval] data)
    X.a = X.raw(1:size(X.raw,1)/2,:);
    X.b = X.raw((size(X.raw,1)/2)+1:end,:);
    
    % get stimulus values for encoding-visual
    stim_vals = stim_details.stimulus(stim_details.encoding==1&stim_details.modality==1);
    
    % split into A and B
    X.sav = stim_vals(1:48);
    X.sbv = stim_vals(49:96);
    
    % get ordered index
    [~,X.sai] = sort(X.sav);
    [~,X.sbi] = sort(X.sbv);
    
    % re-organise regressors
    X.a(:,1:48) = X.a(:,X.sai);
    X.b(:,49:96) = X.b(:,X.sbi+48);
    
    % --- prepare patterns --- %
    % get activation matrix (Y) and 
    Y.raw = patterns(scan_details.encoding==1&scan_details.modality==1,:);
    
    % split into two groups (Ya and Yb)
    Y.a = Y.raw(1:size(Y.raw,1)/2,:);
    Y.b = Y.raw((size(Y.raw,1)/2)+1:end,:);
    
    % --- remove singular dimensions --- %
    % find zero-value rows
    X.a_badRow = all(X.a(:,1:96)==0,2);
    X.b_badRow = all(X.b(:,1:96)==0,2);
    
    % remove zero-value rows
    X.a(X.a_badRow,:) = [];
    X.b(X.b_badRow,:) = [];
    
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
    save([dir_subj,'/rsa-percept-audio/',subj_handle,'_task-percept_rsa-FormattedVolume.mat'],'X','Y')    
    clear X Y i j fn patterns subj_handle dir_subj
    
    fprintf('Subject %02.0f of %02.0f prepared...\n',subj,n_subj)
end

%% Define Searchlight Characteristics
% get the radius of searchlight in voxels
voxRadInSearchLight = scan_search ./ scan_vox;

% define distance from searchlight centre to perimeter in voxels
dist2Perimeter = floor(voxRadInSearchLight);

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
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-percept-audio/',subj_handle,'_task-percept_rsa-FormattedVolume.mat'],'X','Y')
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-percept-audio/',subj_handle,'_task-rf_rsa-mask.mat'],'mask_idx')
    
    % predefine model rdm
    model_rdm{1} = zeros(48,48)+1;
    
    % set items that belong to the same category to -1
    model_rdm{1}(1:12,1:12)   = -1;
    model_rdm{1}(13:24,13:24) = -1;
    model_rdm{1}(25:36,25:36) = -1;
    model_rdm{1}(37:48,37:48) = -1;
    
    % set diagonal to zero (so it can be collapsed)
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
                        voxIdx(end+1) = sub2ind(scan_fov,tmpX(xi),tmpY(yi),tmpZ(zi)); %#ok<AGROW>
                    end
                end
            end
        end
        
        % remove those indices which are not included in mask
        voxIdx(~ismember(voxIdx,M)) = []; %#ok<AGROW>
               
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
    RDM_ldt = rsa.stat.fisherDiscrTRDM_searchlight(X.a,Y.a,X.b,Y.b,1:48,sl_vox,model_rdm);

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
    V.fileprefix = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa-percept-audio/',subj_handle,'_task-percept_rsa-searchlightVisual'];
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
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample    = {[dir_root,'bids_data/derivatives/',subj_handle,'/rsa-percept-audio/',subj_handle,'_task-percept_rsa-searchlightVisual.nii,1']};

    % smooth
    matlabbatch{2}.spm.spatial.smooth.fwhm                      = [8 8 8];
    matlabbatch{2}.spm.spatial.smooth.dtype                     = 0;
    matlabbatch{2}.spm.spatial.smooth.im                        = 0;
    matlabbatch{2}.spm.spatial.smooth.prefix                    = 's';
    matlabbatch{2}.spm.spatial.smooth.data                      = {[dir_root,'bids_data/derivatives/',subj_handle,'/rsa-percept-audio/w',subj_handle,'_task-percept_rsa-searchlightVisual.nii,1']};
    
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
    rMapFiles{subj,1}  = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa-percept-audio/sw',subj_handle,'_task-percept_rsa-searchlightVisual.nii,1'];
end

% make group directory
mkdir([dir_root,'bids_data/derivatives/group/rsa-percept-audio/'])

% delete spm if it exists
if exist([dir_root,'bids_data/derivatives/group/rsa-percept-audio/SPM.mat'],'file'); delete([dir_root,'bids_data/derivatives/group/rsa-percept-audio/SPM.mat']); end

% create second-level glms
matlabbatch{1}.spm.stats.factorial_design.dir                       = {[dir_root,'bids_data/derivatives/group/rsa-percept-audio/']};
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
matlabbatch{1}.spm.stats.fmri_est.spmmat            = {[dir_root,'bids_data/derivatives/group/rsa-percept-audio/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals   = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1;

spm_jobman('run',matlabbatch)
clear matlabbatch

% define contrasts
matlabbatch{1}.spm.stats.con.spmmat(1)                  = {[dir_root,'bids_data/derivatives/group/rsa-percept-audio/SPM.mat']};   
matlabbatch{1}.spm.stats.con.delete                     = 0;    
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name       = 'within>between';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec     = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none';  

spm_jobman('run',matlabbatch)
clear matlabbatch

%% Extract Metrics for Visualisation
% combine visual cluster and save
combine_spm_cluster([dir_root,'bids_data/derivatives/group/rsa-percept-audio/'])

% load SPM details
load([dir_root,'bids_data/derivatives/group/rsa-percept-audio/SPM.mat'],'SPM')

% extract subject values for each cluster
[betas,d] = extract_sample_points([dir_root,'bids_data/derivatives/group/rsa-percept-audio/'],SPM);

% save betas as table
tbl = array2table(betas','VariableNames',{'LeftTemp','RightTemp'});
writetable(tbl,[dir_repos,'data/percept_betas-audio.csv'],'Delimiter',',')

% save effect size as table
tbl = array2table(d','VariableNames',{'LeftTemp','RightTemp'});
writetable(tbl,[dir_repos,'data/percept_cohensD-audio.csv'],'Delimiter',',')

%% --------------------------------------------------------------------- %%
%                                                                         %
%       RESTART FOR RETRIEVAL-BASED AUDITORY ANALYSIS                     %
%                                                                         %
%                                                                         %
% ----------------------------------------------------------------------- %

%% Create GLM Model
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/'])
    mkdir([dir_subj,'rsa-ers/'])
    
    % delete old SPM file if it exists
    if exist([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/SPM.mat'],'file'); delete([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers/SPM.mat']); end
    
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
                case 'TRUMPET';     idx = 1 + offset;
                case 'PIANO';       idx = 2 + offset;
                case 'GUITAR';      idx = 3 + offset;
                case 'ACCORDIAN';   idx = 4 + offset;
                otherwise;          continue
            end
            
            % add key values
            events_onset(count(1,idx),idx) = tbl.onset(e);
            
            % move to select screen
            if strcmpi(tbl.operation{e},'retrieval')
                events_onset(count(1,idx),idx) = events_onset(count(1,idx),idx) + (3 * EEG_sample);
            end
            
            % sort memory
            if tbl.recalled(e) ~= 1 && strcmpi(tbl.operation{e},'retrieval')
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
    else; save([dir_subj,'/rsa-ers-auditory/R.mat'],'R')   
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
        matlabbatch{1}.spm.stats.fmri_spec.dir              = {[dir_subj,'/rsa-ers-auditory']};
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {[dir_subj,'/rsa-ers-auditory/R.mat']}; 
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
    else; matlabbatch{2}.spm.stats.fmri_est.spmmat = {[dir_subj,'/rsa-ers-auditory/SPM.mat']};
    end
    
    % run batch
    spm_jobman('run',matlabbatch)
    clear matlabbatch all_scans bad_scans subj_handle button_onset events_onset trl    
        
    % copy and delete local files if they exist
    if copylocal
        copyfile('C:/tmp_mri/spm/R.mat',[dir_subj,'/rsa-ers-auditory/R.mat'])
        copyfile('C:/tmp_mri/spm/SPM.mat',[dir_subj,'/rsa-ers-auditory/SPM.mat'])
        copyfile('C:/tmp_mri/spm/mask.nii',[dir_subj,'/rsa-ers-auditory/mask.nii'])
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
    nii = load_untouch_nii([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/mask.nii']);
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
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/'])
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/',subj_handle,'_task-rf_rsa-fullVolume.mat'],'scanVec')
    
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
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/',subj_handle,'_task-rf_rsa-maskedVolume.mat'],'patterns')
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/',subj_handle,'_task-rf_rsa-mask.mat'],'mask_idx')
    
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
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/',subj_handle,'_task-rf_rsa-maskedVolume.mat'])
    
    % get mean pattern across all trials and replicate matrix
    meanPattern = repmat(mean(patterns,1),[size(patterns,1), 1]);

    % subtract mean pattern from data
    patterns = patterns - meanPattern;
    
    % save patterns
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/',subj_handle,'_task-rf_rsa-maskedDemeanedVolume.mat'],'patterns')
    
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
    load([dir_subj,'/rsa-ers-auditory/',subj_handle,'_task-rf_rsa-maskedDemeanedVolume.mat'])
    
    % load SPM.mat
    load([dir_subj,'/rsa-ers-auditory/SPM.mat'])
    
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
        block_visual    = double(any(ismember(tbl.modality,'Auditory')));
        
        % cycle through every event
        for e = 1 : size(tbl,1)
            
            % if an event
            if strcmpi(tbl.trial_type(e),'Stimulus Onset')
            
                % add key values
                stim_details.encoding(stim_count) = strcmpi(tbl.operation(e),'encoding');
                stim_details.modality(stim_count) = strcmpi(tbl.modality(e),'auditory');
                
                % get stimulus values
                switch tbl.stimulus{e}
                    case 'TRUMPET';   stim_details.stimulus(stim_count,1) = 1;
                    case 'PIANO';  stim_details.stimulus(stim_count,1) = 2;
                    case 'GUITAR';        stim_details.stimulus(stim_count,1) = 3;
                    case 'ACCORDIAN';        stim_details.stimulus(stim_count,1) = 4;
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
    X.b(X.b_badRow,:) = [];
    
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
    save([dir_subj,'rsa-ers-auditory/',subj_handle,'_task-rf_rsa-FormattedVolume.mat'],'X','Y')
    
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
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    
    % load data
    load([dir_subj,'rsa-ers-auditory/',subj_handle,'_task-rf_rsa-FormattedVolume.mat'])
    load([dir_subj,'rsa-ers-auditory/',subj_handle,'_task-rf_rsa-mask.mat'])
    
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
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/',subj_handle,'_task-rf_rsa-rawRDM.mat'],'RDM_ldt')
    
    % clean up
    clear X Y sl_vox
    
    % add z-value to rdmBrain
    rdmBrain(M(goodSL)) = mean(cat(2,RDM_ldt.ats,RDM_ldt.bts),2);

    % load template nifti
    filename = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/mask.nii'];        
    V = load_untouch_nii(filename);

    % change filename, datatype, and image
    V.fileprefix = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/',subj_handle,'_task-rf_rsa-searchlight_response'];
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
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample    = {[dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/',subj_handle,'_task-rf_rsa-searchlight_response.nii,1']};

    % smooth
    matlabbatch{2}.spm.spatial.smooth.fwhm                      = [8 8 8];
    matlabbatch{2}.spm.spatial.smooth.dtype                     = 0;
    matlabbatch{2}.spm.spatial.smooth.im                        = 1;
    matlabbatch{2}.spm.spatial.smooth.prefix                    = 's';
    matlabbatch{2}.spm.spatial.smooth.data                      = {[dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/w',subj_handle,'_task-rf_rsa-searchlight_response.nii,1']};
    
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
    rMapFiles{subj,1}  = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa-ers-auditory/sw',subj_handle,'_task-rf_rsa-searchlight_response.nii,1'];
end

% create second-level glms
matlabbatch{1}.spm.stats.factorial_design.dir                       = {[dir_root,'bids_data/derivatives/group/rsa-ers-auditory']};
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
matlabbatch{1}.spm.stats.fmri_est.spmmat            = {[dir_root,'bids_data/derivatives/group/rsa-ers-auditory/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals   = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1;

spm_jobman('run',matlabbatch)
clear matlabbatch

% define contrasts
matlabbatch{1}.spm.stats.con.spmmat(1)                  = {[dir_root,'bids_data/derivatives/group/rsa-ers-auditory/SPM.mat']};   
matlabbatch{1}.spm.stats.con.delete                     = 0;    
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name       = 'within>between';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec     = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none';  

spm_jobman('run',matlabbatch)
clear matlabbatch

%% Extract Metrics for Visualisation
% combine visual cluster and save
combine_spm_cluster([dir_root,'bids_data/derivatives/group/rsa-ers-auditory/'])

% load SPM details
load([dir_root,'bids_data/derivatives/group/rsa-ers-auditory/SPM.mat'],'SPM')

% extract subject values for each cluster
[betas,d] = extract_sample_points([dir_root,'bids_data/derivatives/group/rsa-ers-auditory/'],SPM);

% save betas as table
tbl = array2table(betas','VariableNames',{'LeftFrontal'});
writetable(tbl,[dir_repos,'data/fig1_data/ers_resp_betas.csv'],'Delimiter',',')

% save effect size as table
tbl = array2table(d','VariableNames',{'LeftFrontal'});
writetable(tbl,[dir_repos,'data/fig1_data/ers_resp_cohensD.csv'],'Delimiter',',')


%% --------------------------------------------------------------------- %%
%                                                                         %
%       RUN COMBINED AUDITORY ANALYSIS                                    %
%                                                                         %
%                                                                         %
% ----------------------------------------------------------------------- %

%% --- RUN FMRI ANALYSIS ----------------------------------------------- %%

%% Create GLM Model
% set copylocal
if ispc; copylocal = true;
else; copylocal = false;
end

% cycle through each subject
for subj = 1 : n_subj

    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'derivatives/',subj_handle,'/'];
    mkdir([dir_subj,'rsa-correlation-audio/'])

    % delete old SPM file if it exists
    if exist([dir_subj,'rsa-correlation-audio/SPM.mat'],'file'); delete([dir_subj,'rsa-correlation-audio/SPM.mat']); end

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
        tbl = readtable([dir_root,subj_handle,'/func/',subj_handle,'_task-rf_run-',num2str(run),'_events.tsv'],'FileType','text','Delimiter','\t');

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
            if strcmpi(tbl.modality{e},'auditory')

                % add key values
                events_onset(count,1) = tbl.onset(e);
                count = count + 1;
            end

            % get button press onset (if pressed)
            if ~isnan(tbl.rt(e)); button_onset(end+1,1) = tbl.onset(e) + tbl.rt(e); end          %#ok<AGROW>
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
    else; save([dir_subj,'/rsa-correlation-audio/R.mat'],'R')   
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
        matlabbatch{1}.spm.stats.fmri_spec.dir              = {[dir_subj,'/rsa-correlation-audio']};
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {[dir_subj,'/rsa-correlation-audio/R.mat']}; 
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
    else; matlabbatch{2}.spm.stats.fmri_est.spmmat = {[dir_subj,'/rsa-correlation-audio/SPM.mat']};
    end

    % run batch
    spm_jobman('run',matlabbatch)
    clear matlabbatch all_scans bad_scans subj_handle button_onset events_onset trl    

    % copy and delete local files if they exist
    if copylocal
        copyfile('C:/tmp_mri/spm/R.mat',[dir_subj,'/rsa-correlation-audio/R.mat'])
        copyfile('C:/tmp_mri/spm/SPM.mat',[dir_subj,'/rsa-correlation-audio/SPM.mat'])
        copyfile('C:/tmp_mri/spm/mask.nii',[dir_subj,'/rsa-correlation-audio/mask.nii'])
        rmdir('C:/tmp_mri/','s'); 
    end
end

%% Prepare Masks
% cycle through each subject
for subj = 1 : n_subj

    % define subject name
    subj_handle = sprintf('sub-%02.0f',subj);    

    % prepare deformation batch
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def        = {[dir_root,'derivatives/',subj_handle,'/anat/iy_',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.util.defs.comp{1}.inv.space              = {[dir_root,'',subj_handle,'/anat/',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames             = {[dir_root,'derivatives/group/rsa-percept-audio/grand_cluster_dilated.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.weight             = {''};
    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr    = {[dir_root,'derivatives/',subj_handle,'/masks/']};
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.file           = {[dir_root,'derivatives/',subj_handle,'/func/meanua',subj_handle,'_task-rf_run-1_bold.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.preserve           = 0;
    matlabbatch{1}.spm.util.defs.out{1}.push.fwhm               = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.push.prefix             = '';

    % run
    spm_jobman('run',matlabbatch)
    clear matlabbatch

    % change mask name
    movefile([dir_root,'derivatives/',subj_handle,'/masks/wgrand_cluster_dilated.nii'],...
             [dir_root,'derivatives/',subj_handle,'/masks/rsa-percept.nii'])
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
        filename = [dir_root,'derivatives/',subj_handle,'/func/',...
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
    mkdir([dir_root,'derivatives/',subj_handle,'/rsa-correlation-audio/'])
    save([dir_root,'derivatives/',subj_handle,'/rsa-correlation-audio/',subj_handle,'_task-rf_rsa-fullVolume.mat'],'scanVec')

    % apply mask to scans
    fprintf('Masking data...\n')

    % predefine matrix for mask data
    maskImg = zeros(1, prod(scan_fov));

    % load mask and add to matrix    
    nii_1 = load_untouch_nii([dir_root,'derivatives/',subj_handle,'/rsa-correlation-audio/mask.nii']);
    nii_2 = load_untouch_nii([dir_root,'derivatives/',subj_handle,'/masks/rsa-percept.nii']);
    maskImg(1,:) = reshape(nii_1.img==1&nii_2.img==1,1,[]);

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
    save([dir_root,'derivatives/',subj_handle,'/rsa-correlation-audio/',subj_handle,'_task-all_rsa-maskedVolume.mat'],'patterns')
    save([dir_root,'derivatives/',subj_handle,'/rsa-correlation-audio/',subj_handle,'_task-all_rsa-mask.mat'],'mask_idx')

    % clear excess variables
    clear scanVec nii deadIdx patterns mask_idx subjHandle maskImg
end

%% Mean Pattern Subtraction
% cycle through each subject
for subj = 1 : n_subj

    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % load pattern data
    load([dir_root,'derivatives/',subj_handle,'/rsa-correlation-audio/',subj_handle,'_task-all_rsa-maskedVolume.mat'],'patterns')

    % get mean pattern across all trials and replicate matrix
    meanPattern = repmat(mean(patterns,1),[size(patterns,1), 1]);

    % subtract mean pattern from data
    patterns = patterns - meanPattern;

    % save patterns
    save([dir_root,'derivatives/',subj_handle,'/rsa-correlation-audio/',subj_handle,'_task-all_rsa-maskedDemeanedVolume.mat'],'patterns')

    % clean up
    clear subjHandle meanPattern patterns
end

%% Get Trial BOLD 
% predefine mean bold
mean_bold = nan(n_subj,n_trials);

% cycle through each subject
for subj = 1 : n_subj

    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'derivatives/',subj_handle,'/'];

    % load pattern data
    load([dir_subj,'/rsa-correlation-audio/',subj_handle,'_task-all_rsa-maskedDemeanedVolume.mat'],'patterns')

    % load SPM.mat
    load([dir_subj,'/rsa-correlation-audio/SPM.mat'],'SPM')

    % get design matrix (X) and patterns (Y)
    X = SPM.xX.X;
    Y = patterns;

    % cycle through each trial
    for trl = 1 : n_trials

        % get scans invovled in trial
        idx = X(:,trl) > 0;
        x = repmat(X(idx,trl),[1 size(Y,2)]);
        y = Y(idx,:);

        % get mean bold signal
        mean_bold(subj,trl) = mean(mean(x .* y,2),1);
    end
end

% save
mkdir([dir_root,'derivatives/group/rsa-correlation-audio/'])
save([dir_root,'derivatives/group/rsa-correlation-audio/group_task-all_fmri-meanbold'],'mean_bold')
clear X Y x y

%% Calculate LDt
% cycle through each subject
for subj = 1 : n_subj

    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'derivatives/',subj_handle,'/'];

    % load pattern data
    load([dir_subj,'/rsa-correlation-audio/',subj_handle,'_task-all_rsa-maskedDemeanedVolume.mat'],'patterns')

    % load SPM.mat
    load([dir_subj,'/rsa-correlation-audio/SPM.mat'],'SPM')

    % get stimulus tables
    [stim_details,scan_details] = get_auditory_tables(dir_root,subj_handle);

    % get design matrix (X)
    X     = [];
    X.raw = SPM.xX.X;

    % remove scans/regressors that are not auditory (to
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
    Y = [];
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

    % calculate trial-wise linear discriminant T
    RDM.ldtB = rsa.stat.fisherDiscrTRDM_trainAtestB(X.at,Y.a,X.bt,Y.b,size(X.at,2)-96,1:numel(Y.sa));
    RDM.ldtA = rsa.stat.fisherDiscrTRDM_trainAtestB(X.bt,Y.b,X.at,Y.a,size(X.at,2)-96,1:numel(Y.sb));

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
    save([dir_subj,'/rsa-correlation-audio/',subj_handle,'_task-percept_rsa-rdm.mat'],'RDM')
    clear X Y RDM

    fprintf('Subject %02.0f of %02.0f complete...\n',subj,n_subj)
end

%% Extract Vector of Similarity Indices for Each Trial [FORMATING FROM HERE]
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
        X(:,3) = freq{mask}.powspctrm(:,:,:,1);
        X(:,4) = freq{mask}.powspctrm(:,:,:,2);
        X(:,5) = bold_tmp;
        X(:,6) = freq{mask}.trialinfo(:,2);

        % if retrieval, drop forgotten
        if mask == 2; X(freq{mask}.trialinfo(:,1)==0,:) = []; end
        
        % create table
        tbl = array2table(X,'VariableNames',{'rsa','erd','postpow','prepow','bold','conf'});

        % --- fit ERD model --- %
        % define predictors
        preds = {'erd','bold','conf'};
        
        % run linear model
        B = fitlm(tbl,'ResponseVar','rsa',...
                      'RobustOpts','on',...
                      'PredictorVars',preds);

        % extract predictor results
        for p = 1 : numel(preds); grand_freq{mask,1}.(preds{p})(subj,1) = B.Coefficients.tStat(preds{p}); end   
        
        % --- fit pre/post model --- %
        % define predictors
        preds = {'prepow','postpow','bold','conf'};
        
        % run linear model
        B = fitlm(tbl,'ResponseVar','rsa',...
                      'RobustOpts','on',...
                      'PredictorVars',preds);

        % extract predictor results
        for p = 1 : numel(preds); grand_freq{mask,2}.(preds{p})(subj,1) = B.Coefficients.tStat(preds{p}); end   
        
        % --- fit median split model --- %
        % define predictors
        preds = {'prepow','postpow','bold','conf'};
        
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
tail{1}  = [-1 0 0];
tail{2}  = [0 -1 0 0];
tail{3}  = [0 -1 0 0];

% cycle through each model
for model = 1 : 3

    % predefine stat/tbl outputs
    sout = cell(numel(preds),1);
    tout = cell(numel(preds),1);

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

    % add labels
    opt = repmat({'encoding','retrieval'},[1 numel(preds)])';
    var = repmat(preds',[2 1]);
    tbl = addvars(tbl,var(:),'Before','t');
    tbl = addvars(tbl,opt(:),'Before','t');
    disp(tbl)
end

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
        
        % get slope/intercept
        A = [ones(size(freq{mask}.fractal)) freq{mask}.fractal];
        b = freq{mask}.freq;
        
        
        % create regressor matrix
        X      = zeros(numel(trl_nums),4);
        X(:,1) = rsa_vec(subj,mask,trl_nums);
        X(:,2) = freq{mask}.powspctrm;
        X(:,3) = ;
        X(:,4) = bold_tmp;
        X(:,5) = freq{mask}.trialinfo(:,2);

        % if retrieval, drop forgotten
        if mask == 2; X(freq{mask}.trialinfo(:,1)==0,:) = []; end
        
        % create table
        tbl = array2table(X,'VariableNames',{'rsa','erd','frac','bold','conf'});

        % define predictors
        preds = {'erd','frac','bold','conf'};
        
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
tail  = [-1 -1 0 0];

% predefine stat/tbl outputs
sout = cell(numel(preds),1);
tout = cell(numel(preds),1);

% cycle through each parameter
for i = 1 : numel(preds)
    
    % predefine cell for statistics
    cfg             = [];
    cfg.tail        = [0 0] + tail(i);
    cfg.parameter   = preds{i};
    cfg.statistic   = '';
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

%%


% plot
figure; hold on
for i = 1 : size(grand_freq,2)
    subplot(2,1,i); hold on
    for j = 1 : size(grand_freq,1)
        plot(j,grand_freq{j,i}.powspctrm,'k*')
        plot(j+.1,mean(grand_freq{j,i}.powspctrm),'ko')
    end
    set(gca,'xtick',1:size(grand_freq,1),'xticklabel',{'pre','post','conf','bold','acc'});
    xlim([0.5 size(grand_freq,1)+.5])
end

% save data
save([dir_bids,'derivatives/group/rsa-correlation/group_task-retrieval_comb-stat.mat'],'stat','tbl','stat2','tbl2');

% run Hotelling's T [derived from: https://uk.mathworks.com/matlabcentral/fileexchange/2844-hotellingt2] 
for i = 1 : 2
    
    X = cat(2,r(:,1,i),r(:,2,i));
    pval(i) = multi_ttest(X);
end

disp(pval)

%% Extract Raw Power of Cluster 
% prepare table for stat values
tbl = array2table(zeros(n_subj,6),'VariableNames',{'perception','retrieval','forgotten','per_noBold','ret_noBold','ret_noConf'});

% create table
tbl.perception(:,1) = grand_freq{1}.powspctrm;
tbl.retrieval(:,1)  = grand_freq{2}.powspctrm;
tbl.forgotten(:,1)  = grand_freq{3}.powspctrm;
tbl.per_noBold(:,1) = grand_freq_partial{1}.powspctrm;
tbl.ret_noBold(:,1) = grand_freq_partial{2}.powspctrm;
tbl.ret_noConf(:,1) = grand_freq_partial{3}.powspctrm;

% write table
writetable(tbl,[dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-cluster.csv'],'Delimiter',',')
copyfile([dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-cluster.csv'],[dir_repos,'data/fig3_data/group_task-all_eeg-cluster.csv'])

%% Get Whole Brain Map
% get indices of clusters
%clus_idx = stat{1}.negclusterslabelmat==1;

% create source data structure
source                  = [];
source.inside           = stat{1}.inside;
source.dim              = stat{1}.dim;
source.pos              = stat{1}.pos*10;
source.unit             = 'mm';

% define powspctrm of cluster
source.pow              = nan(size(stat{1}.pos,1),1);     
%source.pow(clus_idx)	= stat{1}.stat(clus_idx); 
source.pow	= stat{5}.stat; 
%source.pow(mask_roi{1}.inside)	= stat{5}.stat(mask_roi{1}.inside); 

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
cfg.filename      = [dir_bids,'derivatives/group/rsa-correlation/group_task-ers_eeg-map.nii'];  % enter the desired file name
cfg.filetype      = 'nifti';
cfg.coordsys      = 'spm';
cfg.vmpversion    = 2;
cfg.datatype      = 'float';
ft_volumewrite(cfg,source);      % be sure to use your interpolated source data

% reslice to 1mm isometric to match template MRI
reslice_nii([dir_bids,'derivatives/group/rsa-correlation/group_task-ers_eeg-map.nii'],[dir_bids,'derivatives/group/rsa-correlation/group_task-ers_eeg-map.nii'],[1 1 1]);
copyfile([dir_bids,'derivatives/group/rsa-correlation/group_task-ers_eeg-map.nii'],[dir_repos,'data/fig3_data/group_task-ers_eeg-map.nii'])    

end


