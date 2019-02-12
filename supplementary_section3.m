%% EEG-fMRI: Univariate Analyses
% initialisation
clearvars
clc

% initalise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% add RSA toolbox to path
addpath(genpath('Y:/projects/general/rsatoolbox-develop'))

%% Define Key Parameters
dir_root    = 'Y:/projects/reinstatement_fidelity/';        % data directory
dir_repos   = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
dir_bids    = 'Y:/projects/reinstatement_fidelity/bids_data/';
dir_tool    = 'Y:/projects/general/';
n_subj      = 21;                                           % number of subjects
n_trials    = 192;                                          % number of trials
n_volumes   = 255;                                          % number of volumes
n_runs      = 8;
TR          = 2;
EEG_sample  = 5000;
scan_fov    = [64 64 32];                                   % scan field of view
scan_vox    = [3 3 4];                                      % scan voxel size

% add subfunctions
addpath([dir_repos,'subfunctions'])

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

% load EEG stat
load([dir_bids,'derivatives/group/eeg/group_task-all_eeg-stat.mat'])

% define channels of interest
coi = stat{2}.negclusterslabelmat(roi.inside)==1;

%% Run First-Level Analysis
% cycle through each subject
for subj = 1 : n_subj
    
    % --- Prepare Subject --- %
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_fmri    = [dir_root,'bids_data/',subj_handle,'/'];
    dir_eeg     = [dir_root,'bids_data/sourcedata/',subj_handle,'/'];
    
    % make directory for SPM data
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/univariateEEGfMRI/'])
    
    % delete old SPM file if it exists
    if exist([dir_root,'bids_data/derivatives/',subj_handle,'/univariateEEGfMRI/SPM.mat'],'file'); delete([dir_root,'bids_data/derivatives/',subj_handle,'/univariateEEGfMRI/SPM.mat']); end
    
    % --- Prepare EEG --- %    
    % load data
    load([dir_eeg,'/eeg/',subj_handle,'_task-rf_eeg-source.mat'])
    
    % find audio/video hit/miss trials
    isvis   = nan(size(source.trialinfo));
    isenc   = nan(size(source.trialinfo));
    ismem   = nan(size(source.trialinfo));
    trl_no  = nan(size(source.trialinfo));
    for trl = 1 : numel(source.trialinfo)
        isvis(trl,1) = strcmpi(source.trialinfo{trl}.modality,'visual');
        isenc(trl,1) = strcmpi(source.trialinfo{trl}.operation,'encoding');
        ismem(trl,1) = source.trialinfo{trl}.recalled;
        trl_no(trl,1) = source.trialinfo{trl}.trl_at_retrieval;
    end
    
    % get time-frequency representation of data
    cfg             = [];
    cfg.trials      = isvis & ~isenc;
    cfg.channel     = source.label(coi);
    cfg.keeptrials  = 'yes';
    cfg.method      = 'wavelet';
    cfg.width       = 6;
    cfg.output      = 'pow';
    cfg.toi         = 0.5:0.1:1.5;
    cfg.foi         = 8:2:30; % 100hz sampling
    cfg.pad         = 'nextpow2';
    freq            = ft_freqanalysis(cfg, source); clear source

    % average over time and frequency
    cfg             = [];
    cfg.avgoverchan = 'yes';
    cfg.avgovertime = 'yes';
    cfg.avgoverfreq = 'yes';
    freq            = ft_selectdata(cfg,freq);
       
    % split into hits and misses
    pow_hit             = zscore(freq.powspctrm(ismem(isvis&~isenc)==1));
    pow_miss            = zscore(freq.powspctrm(ismem(isvis&~isenc)==0));
    num_hit             = trl_no(ismem==1&isvis&~isenc);
    num_miss            = trl_no(ismem==0&isvis&~isenc);
    
    % fix trial numbers
    num_hit(num_hit>48 & num_hit<=96)    = num_hit(num_hit>48 & num_hit<=96) - 48;
    num_hit(num_hit>96 & num_hit<=144)   = num_hit(num_hit>96 & num_hit<=144) - 48;
    num_hit(num_hit>144 & num_hit<=192)  = num_hit(num_hit>144 & num_hit<=192) - 96;
    num_miss(num_miss>48 & num_miss<=96)    = num_miss(num_miss>48 & num_miss<=96) - 48;
    num_miss(num_miss>96 & num_miss<=144)   = num_miss(num_miss>96 & num_miss<=144) - 48;
    num_miss(num_miss>144 & num_miss<=192)  = num_miss(num_miss>144 & num_miss<=192) - 96;
    
    % get glm verison of power
    glm_pow             = zeros(96,2);
    glm_pow(num_hit,1)  = pow_hit;
    glm_pow(num_miss,2) = pow_miss;
    
    % --- Prepare fMRI --- %
    % create cell to hold events data
    events_table  = array2table(zeros(n_trials*2,4), 'VariableNames', {'onset','isVisual','isEncoding','isRemembered'});
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
            events_table.onset(count)           = tbl.onset(e);
            events_table.isVisual(count)        = strcmpi(tbl.modality(e),'visual');
            events_table.isEncoding(count)      = strcmpi(tbl.operation(e),'encoding');
            events_table.isRemembered(count)    = tbl.recalled(e);
            
            % get button press onset (if pressed)
            if ~isnan(tbl.rt(e)); button_onset(end+1,1) = tbl.onset(e) + tbl.rt(e); end %#ok<SAGROW>
                
            % update counter
            count = count + 1;
            
        end
    end
       
    % convert event onsets and button presses from EEG samples to seconds
    events_table.onset =  events_table.onset ./ EEG_sample;
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
    mkdir('C:/tmp_mri/spm/'); 
    save('C:/tmp_mri/spm/R.mat','R')
        
    % get all scans for GLM
    all_scans = get_functional_files([dir_root,'bids_data/derivatives/',subj_handle,'/'],'swua',true);
    all_scans = all_scans{1};
    
    % remove bad scans
    all_scans(bad_scans) = [];
    
    % define parameters for GLM
    matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
    matlabbatch{1}.spm.stats.fmri_spec.dir              = {'C:/tmp_mri/spm'};
    matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi       = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress     = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf         = 128;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans       = all_scans;    
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {'C:/tmp_mri/spm/R.mat'};    
    matlabbatch{1}.spm.stats.fmri_spec.timing.units     = 'secs';
    matlabbatch{1}.spm.stats.fmri_spec.timing.RT        = 2;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t    = 32;
    matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0   = 16;     
                
    % add details for each condition
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).name        = 'Visual_Encoding_Recalled';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).onset       = events_table.onset(events_table.isVisual == 1 & events_table.isEncoding == 1 & events_table.isRemembered == 1);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod        = struct('name',{},'param',{},'poly',{});
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name        = 'Visual_Encoding_Forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset       = events_table.onset(events_table.isVisual == 1 & events_table.isEncoding == 1 & events_table.isRemembered == 0);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod        = struct('name',{},'param',{},'poly',{});
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name        = 'Visual_Retrieval_Recalled';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset       = events_table.onset(events_table.isVisual == 1 & events_table.isEncoding == 0 & events_table.isRemembered == 1);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod        = struct('name','pow','param',glm_pow(events_table.isRemembered(events_table.isVisual == 1&events_table.isEncoding==0)==1,1),'poly',1);
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name        = 'Visual_Retrieval_Forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset       = events_table.onset(events_table.isVisual == 1 & events_table.isEncoding == 0 & events_table.isRemembered == 0);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod        = struct('name','pow','param',glm_pow(events_table.isRemembered(events_table.isVisual == 1&events_table.isEncoding==0)==0,2),'poly',1);
       
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name        = 'Auditory_Encoding_Recalled';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset       = events_table.onset(events_table.isVisual == 0 & events_table.isEncoding == 1 & events_table.isRemembered == 1);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod        = struct('name',{},'param',{},'poly',{});
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name        = 'Auditory_Encoding_Forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset       = events_table.onset(events_table.isVisual == 0 & events_table.isEncoding == 1 & events_table.isRemembered == 0);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod        = struct('name',{},'param',{},'poly',{});
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).name        = 'Auditory_Retrieval_Recalled';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).onset       = events_table.onset(events_table.isVisual == 0 & events_table.isEncoding == 0 & events_table.isRemembered == 1);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).pmod        = struct('name',{},'param',{},'poly',{});
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).name        = 'Auditory_Retrieval_Forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).onset       = events_table.onset(events_table.isVisual == 0 & events_table.isEncoding == 0 & events_table.isRemembered == 0);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).pmod        = struct('name',{},'param',{},'poly',{});
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).name        = 'ButtonPress';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).onset       = button_onset;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).duration    = 0.5;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).pmod        = struct('name',{},'param',{},'poly',{});
       
    % --- do parametric for 3 and 4
    visret_onset = events_table.onset(events_table.isVisual == 1 & events_table.isEncoding == 0 & events_table.isRemembered == 1);
    
    % specify model
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
    % estimate GLM
    matlabbatch{1}.spm.stats.fmri_est.write_residuals           = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical          = 1;
    matlabbatch{1}.spm.stats.fmri_est.spmmat                    = {'C:/tmp_mri/spm/SPM.mat'};
      
    % contrast betas
    matlabbatch{2}.spm.stats.con.delete                         = 0;
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.name           = 'retrieval_hits';
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.convec         = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.name           = 'retrieval_misses';
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.convec         = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.spmmat(1)                      = {'C:/tmp_mri/spm/SPM.mat'};    
    
    % run job
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
    % copy data   
    copyfile('C:/tmp_mri/spm/',[dir_root,'bids_data/derivatives/',subj_handle,'/univariateEEGfMRI/'])
    rmdir('C:/tmp_mri/','s');         
end

%% Run Second-Level Stats
% define each first-level contrast name
contrast_labels = {'retrieval_hits' 'retrieval_misses'};

% cycle through each first-level contrast
for i = 1 : numel(contrast_labels)
    
    % define directory for second-level analysis
    stat_dir = [dir_root,'bids_data/derivatives/group/univariateEEGfMRI/',contrast_labels{i},'/'];
    mkdir(stat_dir)
    
    % predefine cell to hold contrast filenames
    contrast_files = cell(n_subj,1);
    
    % cycle through each subkect
    for subj = 1 : n_subj
        
        % add filename of contrast image to group cell
        subj_handle = sprintf('sub-%02.0f',subj);
        contrast_files{subj,1} = [dir_root,'bids_data/derivatives/',subj_handle,'/univariateEEGfMRI/con_',sprintf('%04.0f',i),'.nii'];
        
        clear subj_handle
    end
    
    % specify model
    matlabbatch{1}.spm.stats.factorial_design.cov                       = struct('c',{},'cname',{},'iCFI',{},'iCC',{});
    matlabbatch{1}.spm.stats.factorial_design.multi_cov                 = struct('files',{},'iCFI',{},'iCC',{});
    matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none        = 1;
    matlabbatch{1}.spm.stats.factorial_design.masking.im                = 0;
    matlabbatch{1}.spm.stats.factorial_design.masking.em                = {''};
    matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;
    matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm           = 1;
    matlabbatch{1}.spm.stats.factorial_design.dir                       = {stat_dir};
    matlabbatch{1}.spm.stats.factorial_design.des.t1.scans              = contrast_files;
    
    % estimate model
    matlabbatch{2}.spm.stats.fmri_est.write_residuals                   = 0;
    matlabbatch{2}.spm.stats.fmri_est.method.Classical                  = 1;
    matlabbatch{2}.spm.stats.fmri_est.spmmat                            = {[stat_dir,'SPM.mat']};

    % define contrasts
    matlabbatch{3}.spm.stats.con.delete                                 = 0;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.name                   = 'positive';
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec                 = 1;
    matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep                = 'none';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.name                   = 'negative';
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec                 = -1;
    matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep                = 'none';    
    matlabbatch{3}.spm.stats.con.spmmat(1)                              = {[stat_dir,'SPM.mat']};  

    spm_jobman('run',matlabbatch)
    clear matlabbatch            
end

%% Extract Metrics for Visualisation
% define cluster regions
clus_name{1} = {'occipital','leftTemporalPole','rightTemporalPole'};
clus_name{2} = {'rightFusiform','leftFusiform'};
clus_name{3} = {'occipital','limbic'};

% cycle through each first-level contrast
for i = 1 : 3
    
    % combine visual cluster and save
    combine_spm_cluster([dir_root,'bids_data/derivatives/group/spm/',contrast_labels{i},'/'])

    % load SPM details
    load([dir_root,'bids_data/derivatives/group/spm/',contrast_labels{i},'/SPM.mat'])

    % extract subject values for each cluster
    [betas,d] = extract_sample_points([dir_root,'bids_data/derivatives/group/spm/',contrast_labels{i},'/'],SPM);

    % save betas as table
    tbl = array2table(betas','VariableNames',clus_name{i});
    writetable(tbl,[dir_repos,'data/sup1_data/',contrast_labels{i},'_betas.csv'],'Delimiter',',')

    % save effect size as table
    tbl = array2table(d','VariableNames',clus_name{i});
    writetable(tbl,[dir_repos,'data/sup1_data/',contrast_labels{i},'_cohensD.csv'],'Delimiter',',')
end