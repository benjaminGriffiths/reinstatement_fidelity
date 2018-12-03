%% fMRI: Univariate Analyses
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

% add subfunctions
addpath([dir_repos,'subfunctions'])

%% Run First-Level Analysis
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/',subj_handle,'/'];
    
    % make directory for SPM data
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/spm/'])
    
    % load movement parameters
    R = load([dir_root,'bids_data/derivatives/',subj_handle,'/func/rp_a',subj_handle,'_task-rf_run-1_bold.txt']);
    
    % add regressors that model the per-block linear change
    R(1+(n_volumes*0):n_volumes*1,7) = linspace(0,1,n_volumes);
    R(1+(n_volumes*2):n_volumes*3,8) = linspace(0,1,n_volumes);
    R(1+(n_volumes*4):n_volumes*5,9) = linspace(0,1,n_volumes);
    R(1+(n_volumes*6):n_volumes*7,10) = linspace(0,1,n_volumes);
    R(1+(n_volumes*1):n_volumes*2,11) = linspace(0,1,n_volumes);
    R(1+(n_volumes*3):n_volumes*4,12) = linspace(0,1,n_volumes);
    R(1+(n_volumes*5):n_volumes*6,13) = linspace(0,1,n_volumes);
    R(1+(n_volumes*7):n_volumes*8,14) = linspace(0,1,n_volumes);
    
    % add regressors that model the per-block constant
    R(1+(n_volumes*0):n_volumes*1,15) = zeros(1,n_volumes);
    R(1+(n_volumes*2):n_volumes*3,16) = zeros(1,n_volumes);
    R(1+(n_volumes*4):n_volumes*5,17) = zeros(1,n_volumes);
    R(1+(n_volumes*6):n_volumes*7,18) = zeros(1,n_volumes);
    R(1+(n_volumes*1):n_volumes*2,19) = zeros(1,n_volumes);
    R(1+(n_volumes*3):n_volumes*4,20) = zeros(1,n_volumes);
    R(1+(n_volumes*5):n_volumes*6,21) = zeros(1,n_volumes);
    R(1+(n_volumes*7):n_volumes*8,21) = zeros(1,n_volumes); 
    
    % save nuisance regressors
    save([dir_root,'bids_data/derivatives/',subj_handle,'/spm/R.mat'],'R')        
    
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
        tbl = readtable([dir_root,'bids_data/',subj_handle,'/func/',subj_handle,'_task-rf_run-1_events.tsv'],'FileType','text','Delimiter','\t');
        
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
       
    % define parameters for GLM
    matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
    matlabbatch{1}.spm.stats.fmri_spec.dir              = {[dir_root,'bids_data/derivatives/',subj_handle,'/spm']};
    matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi       = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress     = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf         = 128;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans       = get_functional_files([dir_root,'bids_data/derivatives/',subj_handle,'/'],'swua');
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans       = matlabbatch{1}.spm.stats.fmri_spec.sess.scans{1};
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {[dir_root,'bids_data/derivatives/',subj_handle,'/spm/R.mat']};   
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
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(3).pmod        = struct('name',{},'param',{},'poly',{});
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).name        = 'Visual_Retrieval_Forgotten';
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).onset       = events_table.onset(events_table.isVisual == 1 & events_table.isEncoding == 0 & events_table.isRemembered == 0);
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).duration    = 3;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).tmod        = 0;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).orth        = 1;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(4).pmod        = struct('name',{},'param',{},'poly',{});
    
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
        
    % specify model
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
    % estimate GLM
    matlabbatch{1}.spm.stats.fmri_est.write_residuals           = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical          = 1;
    matlabbatch{1}.spm.stats.fmri_est.spmmat                    = {[dir_root,'bids_data/derivatives/',subj_handle,'/spm/SPM.mat']};
      
    % contrast betas
    matlabbatch{2}.spm.stats.con.delete                         = 0;
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.name           = 'contrastModality_atEncoding';
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.convec         = [1 1 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.name           = 'contrastModality_atRetrieval';
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.convec         = [0 0 1 1 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{4}.tcon.name           = 'contrastModality_retrievalSuccess';
    matlabbatch{2}.spm.stats.con.consess{4}.tcon.convec         = [0 0 1 0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{4}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{4}.tcon.name           = 'retrievalSuccess_visual';
    matlabbatch{2}.spm.stats.con.consess{4}.tcon.convec         = [0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{4}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{5}.tcon.name           = 'retrievalSuccess_auditory';
    matlabbatch{2}.spm.stats.con.consess{5}.tcon.convec         = [0 0 0 0 0 0 1 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{5}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.spmmat(1)                      = {[dir_root,'bids_data/derivatives/',subj_handle,'/spm/SPM.mat']};    
    
    % run job
    spm_jobman('run',matlabbatch)
    clear matlabbatch
        
end