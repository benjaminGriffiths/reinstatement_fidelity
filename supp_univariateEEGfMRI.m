function supp_univariateEEGfMRI

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
dir_tool    = 'Y:/projects/general/';
n_subj      = 21;                                           % number of subjects
n_trials    = 192;                                          % number of trials
n_volumes   = 255;                                          % number of volumes
n_runs      = 8;
TR          = 2;
EEG_sample  = 5000;

% add subfunctions
addpath([dir_repos,'subfunctions'])

%% Prepare ROI
% load sourcemodel
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat'],'sourcemodel');

% load whole brain
mri = ft_read_mri([dir_root,'sourcedata/masks/whole_brain.nii']);

% interpolate clusters with grid
cfg             = [];
cfg.parameter	= 'anatomy';
roi             = ft_sourceinterpolate(cfg,mri,sourcemodel);

% define additional roi parameters
roi.inside      = roi.anatomy(:) > 0;
roi.pos         = sourcemodel.pos;

% load EEG stat
load([dir_root,'derivatives/group/eeg/group_task-rf_pow-wavestat.mat'],'stat')

%% Run First-Level Analysis
% cycle through each subject
for subj = 1 : n_subj
    
    % --- Prepare Subject --- %
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_eeg     = [dir_root,'sourcedata/',subj_handle,'/'];
    
    % make directory for SPM data
    mkdir([dir_root,'derivatives/',subj_handle,'/univariateEEGfMRI/'])
    
    % delete old SPM file if it exists
    if exist([dir_root,'derivatives/',subj_handle,'/univariateEEGfMRI/SPM.mat'],'file'); delete([dir_root,'derivatives/',subj_handle,'/univariateEEGfMRI/SPM.mat']); end
    
    % --- Prepare EEG --- %    
    % load data
    load([dir_eeg,'/eeg/',subj_handle,'_task-rf_eeg-source.mat'],'source')
    
    % find audio/video hit/miss trials
    isvis       = nan(size(source.trialinfo));
    isenc       = nan(size(source.trialinfo));
    ismem       = nan(size(source.trialinfo));
    trl_no_enc  = nan(size(source.trialinfo));
    trl_no_ret  = nan(size(source.trialinfo));
    for trl = 1 : numel(source.trialinfo)
        isvis(trl,1) = strcmpi(source.trialinfo{trl}.modality,'visual');
        isenc(trl,1) = strcmpi(source.trialinfo{trl}.operation,'encoding');
        ismem(trl,1) = source.trialinfo{trl}.recalled;
        
        if isenc(trl,1) == 1
            trl_no_enc(trl,1) = source.trialinfo{trl}.trl_at_encoding;
        else
            trl_no_ret(trl,1) = source.trialinfo{trl}.trl_at_retrieval;
        end
    end
    
    % fix encoding trial numbers
    trl_no_enc(trl_no_enc>144 & trl_no_enc<=192)  = trl_no_enc(trl_no_enc>144 & trl_no_enc<=192) + 144;
    trl_no_enc(trl_no_enc>96 & trl_no_enc<=144)   = trl_no_enc(trl_no_enc>96 & trl_no_enc<=144) + 96;
    trl_no_enc(trl_no_enc>48 & trl_no_enc<=96)    = trl_no_enc(trl_no_enc>48 & trl_no_enc<=96) + 48;
    
    % fix retrieval trial numbers
    trl_no_ret(trl_no_ret>144 & trl_no_ret<=192)  = trl_no_ret(trl_no_ret>144 & trl_no_ret<=192) + 192;
    trl_no_ret(trl_no_ret>96 & trl_no_ret<=144)   = trl_no_ret(trl_no_ret>96 & trl_no_ret<=144) + 144;
    trl_no_ret(trl_no_ret>48 & trl_no_ret<=96)    = trl_no_ret(trl_no_ret>48 & trl_no_ret<=96) + 96;
    trl_no_ret(trl_no_ret>0 & trl_no_ret<=48)     = trl_no_ret(trl_no_ret>0 & trl_no_ret<=48) + 48;
    
    % conjoin encoding and retrieval trial numbers
    trl_no = cat(1,trl_no_enc(~isnan(trl_no_enc)),trl_no_ret(~isnan(trl_no_ret)));
    
    % get time-frequency representation of data
    cfg             = [];
    cfg.keeptrials  = 'yes';
    cfg.method      = 'wavelet';
    cfg.width       = 6;
    cfg.output      = 'pow';
    cfg.toi         = [-1:0.1:-0.375, 0.5:0.1:1.5];
    cfg.foi         = 8:2:30; % 100hz sampling
    cfg.pad         = 'nextpow2';
    freq            = ft_freqanalysis(cfg, source); clear source

    % baseline correct
    pre_toi         = freq.time < 0;
    post_toi        = freq.time > 0;
    freq.powspctrm  = mean(freq.powspctrm(:,:,:,post_toi),4) - mean(freq.powspctrm(:,:,:,pre_toi),4); % get post > pre diff
    
    % restrict to relevant channels
    pow = zeros(size(freq.powspctrm,1),1);
    for trl = 1 : size(freq.powspctrm,1)
        
        % define coi based on conditon
        if isvis(trl)
            coi = stat{1}.negclusterslabelmat(stat{1}.inside==1) == 1;
        else
            coi = stat{4}.negclusterslabelmat(stat{1}.inside==1) == 1;
        end
        
        pow(trl,1) = mean(mean(freq.powspctrm(trl,coi,:),3),2);
    end
    freq.powspctrm = pow;
    
    % z-score within condition
    for a = 0:1
        for b = 0:1
            for c = 0:1
                freq.powspctrm(isvis==a & isenc==b & ismem==c) = zscore(freq.powspctrm(isvis==a & isenc==b & ismem==c));
            end
        end
    end
     
    % create glm pow structure
    glm_pow = zeros(n_trials*2,1);
    for i = 1 : numel(freq.powspctrm)
        glm_pow(trl_no(i)) = freq.powspctrm(i);
    end
    
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
            
            % add key values
            events_table.onset(count)           = tbl.onset(e);
            events_table.isVisual(count)        = strcmpi(tbl.modality(e),'visual');
            events_table.isEncoding(count)      = strcmpi(tbl.operation(e),'encoding');
            events_table.isRemembered(count)    = tbl.recalled(e);
            
            % get button press onset (if pressed)
            if ~isnan(tbl.rt(e)); button_onset(end+1,1) = tbl.onset(e) + tbl.rt(e); end %#ok<AGROW>
                
            % update counter
            count = count + 1;
            
        end
    end
       
    % convert event onsets and button presses from EEG samples to seconds
    events_table.onset =  events_table.onset ./ EEG_sample;
    button_onset = button_onset ./ EEG_sample;    
    
    % load movement parameters
    R = load([dir_root,'derivatives/',subj_handle,'/func/rp_a',subj_handle,'_task-rf_run-1_bold.txt']);
    
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
    all_scans = get_functional_files([dir_root,'derivatives/',subj_handle,'/'],'swua',true);
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
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(1).pmod        = struct('name','pow','param',glm_pow(events_table.isVisual == 1 & events_table.isEncoding == 1 & events_table.isRemembered == 1),'poly',1);
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).name        = 'Visual_Encoding_Forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).onset       = events_table.onset(events_table.isVisual == 1 & events_table.isEncoding == 1 & events_table.isRemembered == 0);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(2).pmod        = struct('name','pow','param',glm_pow(events_table.isVisual == 1 & events_table.isEncoding == 1 & events_table.isRemembered == 0),'poly',1);
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).name        = 'Visual_Retrieval_Recalled';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).onset       = events_table.onset(events_table.isVisual == 1 & events_table.isEncoding == 0 & events_table.isRemembered == 1);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod        = struct('name','pow','param',glm_pow(events_table.isVisual == 1 & events_table.isEncoding == 0 & events_table.isRemembered == 1),'poly',1);
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name        = 'Visual_Retrieval_Forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset       = events_table.onset(events_table.isVisual == 1 & events_table.isEncoding == 0 & events_table.isRemembered == 0);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod        = struct('name','pow','param',glm_pow(events_table.isVisual == 1 & events_table.isEncoding == 0 & events_table.isRemembered == 0),'poly',1);
       
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).name        = 'Auditory_Encoding_Recalled';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).onset       = events_table.onset(events_table.isVisual == 0 & events_table.isEncoding == 1 & events_table.isRemembered == 1);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(5).pmod        = struct('name','pow','param',glm_pow(events_table.isVisual == 0 & events_table.isEncoding == 1 & events_table.isRemembered == 1),'poly',1);
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).name        = 'Auditory_Encoding_Forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).onset       = events_table.onset(events_table.isVisual == 0 & events_table.isEncoding == 1 & events_table.isRemembered == 0);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(6).pmod        = struct('name','pow','param',glm_pow(events_table.isVisual == 0 & events_table.isEncoding == 1 & events_table.isRemembered == 0),'poly',1);
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).name        = 'Auditory_Retrieval_Recalled';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).onset       = events_table.onset(events_table.isVisual == 0 & events_table.isEncoding == 0 & events_table.isRemembered == 1);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(7).pmod        = struct('name','pow','param',glm_pow(events_table.isVisual == 0 & events_table.isEncoding == 0 & events_table.isRemembered == 1),'poly',1);
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).name        = 'Auditory_Retrieval_Forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).onset       = events_table.onset(events_table.isVisual == 0 & events_table.isEncoding == 0 & events_table.isRemembered == 0);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(8).pmod        = struct('name','pow','param',glm_pow(events_table.isVisual == 0 & events_table.isEncoding == 0 & events_table.isRemembered == 0),'poly',1);
    
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
    matlabbatch{1}.spm.stats.fmri_est.spmmat                    = {'C:/tmp_mri/spm/SPM.mat'};
      
    % contrast betas
    matlabbatch{2}.spm.stats.con.delete                         = 0;
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.name           = 'retrieval_visual_hits';
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.convec         = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.name           = 'retrieval_audio_hits';
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.convec         = [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{3}.tcon.name           = 'encoding_visual_all';
    matlabbatch{2}.spm.stats.con.consess{3}.tcon.convec         = [1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{3}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{4}.tcon.name           = 'encoding_audio_all';
    matlabbatch{2}.spm.stats.con.consess{4}.tcon.convec         = [0 0 0 0 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{4}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{5}.tcon.name           = 'retrieval_memory';
    matlabbatch{2}.spm.stats.con.consess{5}.tcon.convec         = [0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{5}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{6}.tcon.name           = 'encoding_memory';
    matlabbatch{2}.spm.stats.con.consess{6}.tcon.convec         = [1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{6}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{7}.tcon.name           = 'retrieval_modality';
    matlabbatch{2}.spm.stats.con.consess{7}.tcon.convec         = [0 0 1 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{7}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{8}.tcon.name           = 'encoding_modality';
    matlabbatch{2}.spm.stats.con.consess{8}.tcon.convec         = [1 1 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{8}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.spmmat(1)                      = {'C:/tmp_mri/spm/SPM.mat'};    
    
    % run job
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
    % copy data   
    copyfile('C:/tmp_mri/spm/',[dir_root,'derivatives/',subj_handle,'/univariateEEGfMRI/'])
    rmdir('C:/tmp_mri/','s');         
end

%% Run Second-Level Stats
% define each first-level contrast name
contrast_labels = {'retrieval_visual_hits' 'retrieval_audio_hits' 'encoding_visual_all' 'encoding_audio_all' 'retrieval_memory' 'encoding_memory' 'retrieval_modality' 'encoding_modality'};

% cycle through each first-level contrast
for i = 1 : numel(contrast_labels)
    
    % define directory for second-level analysis
    stat_dir = [dir_root,'derivatives/group/univariateEEGfMRI/',contrast_labels{i},'/'];
    mkdir(stat_dir)
    
    % delete old SPM file if it exists
    if exist([stat_dir,'/SPM.mat'],'file'); delete([stat_dir,'/SPM.mat']); end
    
    % predefine cell to hold contrast filenames
    contrast_files = cell(n_subj,1);
    
    % cycle through each subkect
    for subj = 1 : n_subj
        
        % add filename of contrast image to group cell
        subj_handle = sprintf('sub-%02.0f',subj);
        contrast_files{subj,1} = [dir_root,'derivatives/',subj_handle,'/univariateEEGfMRI/con_',sprintf('%04.0f',i),'.nii'];
        
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
% define each first-level contrast name
% contrast_labels = {'retrieval_memory' 'encoding_memory' 'retrieval_modality' 'encoding_modality'};
% 
% % define cluster regions
% clus_name{1} = {'occipital','parietal'};
% clus_name{2} = {'occipital','parietal'};
% clus_name{3} = {'fusiform','parietal','leftFrontal','rightFrontal','leftIFG'};
% clus_name{4} = {'fusiform','parietal','leftFrontal','rightFrontal'};
% 
% % cycle through each first-level contrast
% for i = 1 : numel(contrast_labels)
%     
%     % combine visual cluster and save
%     combine_spm_cluster([dir_root,'derivatives/group/univariateEEGfMRI/',contrast_labels{i},'/'])
% 
%     % load SPM details
%     load([dir_root,'derivatives/group/univariateEEGfMRI/',contrast_labels{i},'/SPM.mat'],'SPM')
% 
%     % extract subject values for each cluster
%     [betas,d] = extract_sample_points([dir_root,'derivatives/group/univariateEEGfMRI/',contrast_labels{i},'/'],SPM);
% 
%     % save betas as table
%     tbl = array2table(betas','VariableNames',clus_name{i});
%     writetable(tbl,[dir_repos,'data/sup3_data/',contrast_labels{i},'_betas.csv'],'Delimiter',',')
% 
%     % save effect size as table
%     tbl = array2table(d','VariableNames',clus_name{i});
%     writetable(tbl,[dir_repos,'data/sup3_data/',contrast_labels{i},'_cohensD.csv'],'Delimiter',',')
% end