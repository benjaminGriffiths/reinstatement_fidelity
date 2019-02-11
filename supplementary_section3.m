%% EEG-fMRI: Univariate Analyses
% initialisation
clearvars
clc

% initalise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

%% Define Key Parameters
dir_root    = 'Y:/projects/reinstatement_fidelity/';        % data directory
dir_repos   = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
n_subj      = 21;                                           % number of subjects
n_trials    = 192;                                          % number of trials
n_volumes   = 255;                                          % number of volumes
n_runs      = 8;
TR          = 2;
EEG_sample  = 5000;

% add subfunctions
addpath([dir_repos,'subfunctions'])

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
    [audio_bool,hit_bool] = get_splits_eeg(source);
    
    % get time-frequency representation of data
    cfg             = [];
    cfg.keeptrials  = 'yes';
    cfg.method      = 'wavelet';
    cfg.width       = 6;
    cfg.output      = 'pow';
    cfg.toi         = 0.5:0.1:1.5;
    cfg.foi         = 8:2:30; % 100hz sampling
    cfg.pad         = 'nextpow2';
    freq            = ft_freqanalysis(cfg, source); clear source

    % convert powspctrm to single
    freq.powspctrm = single(freq.powspctrm);
    
    % average over time and frequency
    cfg             = [];
    cfg.avgoverchan = 'yes';
    cfg.avgovertime = 'yes';
    cfg.avgoverfreq = 'yes';
    freq            = ft_selectdata(cfg,freq);
       
    % extract power spectrum and z-transform for each condition seperately
    powspctrm = freq.powspctrm;
    powspctrm(hit_bool==1 & audio_bool==1) = zscore(powspctrm(hit_bool==1 & audio_bool==1));
    powspctrm(hit_bool==1 & audio_bool==0) = zscore(powspctrm(hit_bool==1 & audio_bool==0));
    powspctrm(hit_bool==0 & audio_bool==1) = zscore(powspctrm(hit_bool==0 & audio_bool==1));
    powspctrm(hit_bool==0 & audio_bool==0) = zscore(powspctrm(hit_bool==0 & audio_bool==0));
       
    % define a vector the length of all trials (including those rejected
    % due to artifacts)
    glm_pow = zeros(n_trials,1);

    % cycle through each freq trial and add to pow vector
    for i = 1 : numel(freq.powspctrm)
        glm_pow(freq.trialinfo{i}.trlNo_retrieval,1) = powspctrm(i,1);
    end     
    
    clear powspctrm freq source
    
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
    save([dir_root,'bids_data/derivatives/',subj_handle,'/univariateEEGfMRI/R.mat'],'R')            
    
    % get all scans for GLM
    all_scans = get_functional_files([dir_root,'bids_data/derivatives/',subj_handle,'/'],'swua');
    all_scans = all_scans{1};
    
    % remove bad scans
    all_scans(bad_scans) = [];
    
    % define parameters for GLM
    matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
    matlabbatch{1}.spm.stats.fmri_spec.dir              = {[dir_root,'bids_data/derivatives/',subj_handle,'/univariateEEGfMRI']};
    matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi       = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress     = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf         = 128;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans       = all_scans;
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {[dir_root,'bids_data/derivatives/',subj_handle,'/univariateEEGfMRI/R.mat']};   
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
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod        = struct('name','pow','param',glm_pow(events_table.isVisual(events_table.isEncoding==0) == 1 & events_table.isRemembered(events_table.isEncoding==0) == 1),'poly',1);
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name        = 'Visual_Retrieval_Forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset       = events_table.onset(events_table.isVisual == 1 & events_table.isEncoding == 0 & events_table.isRemembered == 0);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod        = struct('name','pow','param',glm_pow(events_table.isVisual(events_table.isEncoding==0) == 1 & events_table.isRemembered(events_table.isEncoding==0) == 0),'poly',1);
    
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
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).pmod        = struct('name','pow','param',glm_pow(events_table.isVisual(events_table.isEncoding==0) == 0 & events_table.isRemembered(events_table.isEncoding==0) == 1),'poly',1);
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).name        = 'Auditory_Retrieval_Forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).onset       = events_table.onset(events_table.isVisual == 0 & events_table.isEncoding == 0 & events_table.isRemembered == 0);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).pmod        = struct('name','pow','param',glm_pow(events_table.isVisual(events_table.isEncoding==0) == 0 & events_table.isRemembered(events_table.isEncoding==0) == 0),'poly',1);
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).name        = 'ButtonPress';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).onset       = button_onset;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).duration    = 0.5;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(9).pmod        = struct('name',{},'param',{},'poly',{});
        
    % specify model
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
    % estimate GLM
    matlabbatch{1}.spm.stats.fmri_est.write_residuals           = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical          = 1;
    matlabbatch{1}.spm.stats.fmri_est.spmmat                    = {[dir_root,'bids_data/derivatives/',subj_handle,'/univariateEEGfMRI/SPM.mat']};
      
    % contrast betas
    matlabbatch{2}.spm.stats.con.delete                         = 0;
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.name           = 'contrastModality_atRetrieval';
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.convec         = [0 0 1 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.name           = 'retrievalSuccess_visual';
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.convec         = [0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{3}.tcon.name           = 'retrievalSuccess_auditory';
    matlabbatch{2}.spm.stats.con.consess{3}.tcon.convec         = [0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{3}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.spmmat(1)                      = {[dir_root,'bids_data/derivatives/',subj_handle,'/univariateEEGfMRI/SPM.mat']};    
    
    % run job
    spm_jobman('run',matlabbatch)
    clear matlabbatch
        
end
