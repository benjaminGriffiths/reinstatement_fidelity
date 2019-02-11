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
    mkdir([dir_subj,'entropy/'])
    
    % delete old SPM file if it exists
    if exist([dir_subj,'entropy/SPM.mat'],'file'); delete([dir_subj,'entropy/SPM.mat']); end
    
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
    else; save([dir_subj,'/entropy/R.mat'],'R')   
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
        matlabbatch{1}.spm.stats.fmri_spec.dir              = {[dir_subj,'/entropy']};
        matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {[dir_subj,'/entropy/R.mat']}; 
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
    else; matlabbatch{2}.spm.stats.fmri_est.spmmat = {[dir_subj,'/entropy/SPM.mat']};
    end
    
    % run batch
    spm_jobman('run',matlabbatch)
    clear matlabbatch all_scans bad_scans subj_handle button_onset events_onset trl    
    
    % copy and delete local files if they exist
    if copylocal
        
        % copy spm files
        copyfile('C:/tmp_mri/spm/R.mat',[dir_subj,'/entropy/R.mat'])
        copyfile('C:/tmp_mri/spm/SPM.mat',[dir_subj,'/entropy/SPM.mat'])
        copyfile('C:/tmp_mri/spm/mask.nii',[dir_subj,'/entropy/mask.nii'])
        
        % copy betas
        files = dir('C:/tmp_mri/spm/beta*');
        for i = 1 : numel(files)
            copyfile(['C:/tmp_mri/spm/',files(i).name],[dir_subj,'entropy/',files(i).name])
        end
        
        % remove directory
        rmdir('C:/tmp_mri/','s'); 
    end
end

%% Calculate Entropy
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    
    % cycle through each mask
    for mask = 1 : numel(mask_names)

        % load mask and add to matrix    
        nii_1 = load_untouch_nii([dir_root,'bids_data/derivatives/',subj_handle,'/entropy/mask.nii']);
        nii_2 = load_untouch_nii([dir_root,'bids_data/derivatives/',subj_handle,'/masks/rsa-',mask_names{mask},'.nii']);
        mask_img = reshape(nii_1.img == 1 & nii_2.img == 1,[],1);
        
        % cycle through each beta
        for i = 1 : n_trials

            % load beta
            nii = load_untouch_nii([dir_subj,'entropy/',sprintf('beta_%04.0f.nii',i)]);
            
            % extract vector of beta
            beta = nii.img(:);
            
            % mask beta
            beta(mask_img~=1) = [];
                       
            % calculate Shannon's entropy
            H(subj,mask,i) = log(-wentropy(beta,'shannon'));            
        end        
    end
end

% save
save([dir_bids,'derivatives/group/entropy/group_task-all_fmri-H.mat'],'H')

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
group_freq   = cell(n_subj,1);

% copy rsa-correlation data if it exists
if exist([dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-freq.mat'],'file')
    copyfile([dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-freq.mat'],[dir_bids,'derivatives/group/entropy/group_task-all_eeg-freq.mat'])
    load([dir_bids,'derivatives/group/entropy/group_task-all_eeg-freq.mat'])
    
% if not, calculate
else
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
    save([dir_bids,'derivatives/group/entropy/group_task-all_eeg-freq.mat'],'group_freq','-v7.3'); 
end

%% Correlate Entropy and Power
% load EEG stat
load([dir_bids,'derivatives/group/eeg/group_task-all_eeg-stat.mat'])

% load in similarity index
load([dir_root,'bids_data/derivatives/group/entropy/group_task-all_fmri-H.mat'])

% cycle through each subject
for subj = 1 : n_subj
    
    % cycle through each mask
    for mask = 1 : numel(mask_names)
    
        % get channels in roi
        coi = stat{mask}.negclusterslabelmat(stat{mask}.inside==1)==1;
        
        % select only remembered items if second mask
        if mask == 2        
            
            % extract recalled trial numbers from powspctrm
            idx         = group_freq{subj,mask}.trialinfo(:,2)==1;
            trl_nums    = group_freq{subj,mask}.trialinfo(idx,1);
            
            % extract relevant entropy values
            Ht = squeeze(H(subj,mask,[1:48 97:144]));
            
        else            
            % extract all trial numbers from powspctrm
            idx         = 1 : size(group_freq{subj,mask}.trialinfo,1);
            trl_nums    = group_freq{subj,mask}.trialinfo(idx,1);
                        
            % extract relevant entropy values
            Ht = squeeze(H(subj,mask,[49:96 145:192]));
        end

        % fix numbers
        trl_nums(trl_nums>48 & trl_nums<=96)    = trl_nums(trl_nums>48 & trl_nums<=96) - 48;
        trl_nums(trl_nums>96 & trl_nums<=144)   = trl_nums(trl_nums>96 & trl_nums<=144) - 48;
        trl_nums(trl_nums>144 & trl_nums<=192)  = trl_nums(trl_nums>144 & trl_nums<=192) - 96;

        % extract vectors for remembered items
        X = Ht(trl_nums);
        Y = nanmean(nanmean(group_freq{subj,mask}.powspctrm(idx,coi,:),3),2);

        % correlate
        r(subj,mask) = atanh(corr(X,Y));
    end
end

% prepare data structure for one-sample t-test
grand_freq{1}              = [];
grand_freq{1}.time         = 1;
grand_freq{1}.freq         = 1;
grand_freq{1}.label        = {'dummy'};
grand_freq{1}.dimord       = 'subj_chan_freq_time';
grand_freq{1}.powspctrm    = r(:,1);

% duplicate data structure for ERS data
grand_freq{2}               = grand_freq{1};
grand_freq{2}.powspctrm     = r(:,2);

% save data
save([dir_bids,'derivatives/group/entropy/group_task-all_comb-freq.mat'],'grand_freq'); 

%% Run Statistics
% predefine cell for statistics
cfg             = [];
cfg.tail        = -1;
[stat,tbl]      = run_oneSampleT(cfg, grand_freq);
   
% save data
save([dir_bids,'derivatives/group/entropy/group_task-retrieval_comb-stat.mat'],'stat','tbl');

%% Extract Raw Power of Cluster 
% prepare table for stat values
tbl = array2table(zeros(n_subj,2),'VariableNames',{'perception','retrieval'});

% create table
tbl.perception(:,1) = grand_freq{1}.powspctrm;
tbl.retrieval(:,1)  = grand_freq{2}.powspctrm;

% write table
writetable(tbl,[dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-cluster.csv'],'Delimiter',',')
copyfile([dir_bids,'derivatives/group/rsa-correlation/group_task-all_eeg-cluster.csv'],[dir_repos,'data/fig3_data/group_task-all_eeg-cluster.csv'])
