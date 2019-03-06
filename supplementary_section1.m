%% fMRI: Univariate Analyses
% initialisation
clearvars
clc

% initalise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

% add RSA toolbox to path
addpath(genpath('Y:/projects/general/rsatoolbox-develop'))
addpath(genpath('Y:\projects\general\GLMdenoise-1.4\'))

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

% add subfunctions
addpath([dir_repos,'subfunctions'])

%% Get Denoised GLM
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/',subj_handle,'/'];
    
    % predefine data and design cells
    data    = cell(1,8);
    design  = cell(1,8);
    
    % cycle through each run
    for run = 1 : n_runs

        % load data
        nii = load_nii([dir_root,'bids_data/derivatives/',subj_handle,'/func/swua',subj_handle,'_task-rf_run-',num2str(run),'_bold.nii']);
        
        % add bold to data cell
        data{run} = single(nii.img(:,:,:,4:end-5));
        
        % load event table
        tbl = readtable([dir_root,'bids_data/',subj_handle,'/func/',subj_handle,'_task-rf_run-',num2str(run),'_events.tsv'],'FileType','text','Delimiter','\t');
        
        % remove first three volumes and last five volumes
        tbl = tbl(4:end-5,:);
        
        % predefine design matrix
        design{run} = zeros(size(data{run},4),6);
        
        % cycle through every event
        for e = 1 : size(tbl,1)
            
            % find stimulus onset
            if strcmpi(tbl.trial_type{e},'stimulus onset')
                
                % get all other trials
                all_others  = tbl.onset(~ismember(1:numel(tbl.onset),e));
                all_labels  = tbl.trial_type(~ismember(1:numel(tbl.onset),e));
                                
                % find nearest scan
                idx         = knnsearch(all_others,tbl.onset(e));
                
                % get volume number
                vol_num = str2double(all_labels{idx}(7:end))-3;
                
                % define regressor num
                if strcmpi(tbl.operation{e},'encoding') && strcmpi(tbl.modality{e},'visual')
                    reg_num = 1;
                elseif strcmpi(tbl.operation{e},'encoding') && strcmpi(tbl.modality{e},'auditory')
                    reg_num = 2;
                elseif strcmpi(tbl.operation{e},'retrieval') && strcmpi(tbl.modality{e},'visual') && tbl.recalled(e) == 1
                    reg_num = 3;
                elseif strcmpi(tbl.operation{e},'retrieval') && strcmpi(tbl.modality{e},'auditory') && tbl.recalled(e) == 1
                    reg_num = 4;
                elseif strcmpi(tbl.operation{e},'retrieval') && strcmpi(tbl.modality{e},'visual') && tbl.recalled(e) == 0
                    reg_num = 5;
                elseif strcmpi(tbl.operation{e},'retrieval') && strcmpi(tbl.modality{e},'auditory') && tbl.recalled(e) == 0
                    reg_num = 6;
                end
                                    
                % add to design
                design{run}(vol_num,reg_num) = 1;    
            end            
        end
    end
           
    % combine runs 1-4 and 5-8
    cat_data{1} = cat(4,data{1},data{2},data{3},data{4});
    cat_data{2} = cat(4,data{5},data{6},data{7},data{8});
    
    % combine designs 1-4 and 5-8
    cat_design{1} = cat(1,design{1},design{2},design{3},design{4});
    cat_design{2} = cat(1,design{5},design{6},design{7},design{8});
    
    % run glm denoise
    [results] = GLMdenoisedata(cat_design,cat_data,3,TR,'assume',[],struct('numboots',100,'numpcstotry',20,'wantparametric',1),[]);
    
    % extract beta maps
    beta = results.parametric.parameters;
    
    % save denoised data
    for run = 1 : n_runs
        
        % get idx
        idx     = mod(run,4); if idx==0;idx=4;end
        block   = floor((run-1)/4)+1;
        trls    = size(design{1})*(idx-1)+1 : size(design{1},1)*idx;
        
        % load old data
        nii = load_nii([dir_root,'bids_data/derivatives/',subj_handle,'/func/swua',subj_handle,'_task-rf_run-',num2str(run),'_bold.nii']);
        
        % replace image with new data
        nii.img = cat_data{block}(:,:,:,trls);
            
        % update template name and variables
        nii.hdr.dime.dim(5) = 247;
        nii.fileprefix = [dir_root,'bids_data/derivatives/',subj_handle,'/func/dswua',subj_handle,'_task-rf_run-',num2str(run),'_bold'];
        
        % save nii
        save_nii(nii,[dir_root,'bids_data/derivatives/',subj_handle,'/func/dswua',subj_handle,'_task-rf_run-',num2str(run),'_bold.nii'])
        
    end
end  
    
%% Run First-Level Analysis
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/',subj_handle,'/'];
    
    % make directory for SPM data
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/spm_denoised/'])
    
    % delete old SPM file if it exists
    if exist([dir_root,'bids_data/derivatives/',subj_handle,'/spm_denoised/SPM.mat'],'file'); delete([dir_root,'bids_data/derivatives/',subj_handle,'/spm_denoised/SPM.mat']); end
    
    % create cell to hold events data
    events_table  = array2table(zeros(n_trials*2,4), 'VariableNames', {'onset','isVisual','isEncoding','isRemembered'});
    count = 1;
    
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
            
            % update counter
            count = count + 1;
            
        end
    end
       
    % convert event onsets and button presses from EEG samples to seconds
    events_table.onset =  events_table.onset ./ EEG_sample;
    
    % get all scans for GLM
    all_scans = get_functional_files([dir_root,'bids_data/derivatives/',subj_handle,'/'],'dswua',true);
    all_scans = all_scans{1};
    
    % define parameters for GLM
    matlabbatch{1}.spm.stats.fmri_spec.volt             = 1;
    matlabbatch{1}.spm.stats.fmri_spec.global           = 'None';
    matlabbatch{1}.spm.stats.fmri_spec.mthresh          = 0.8;
    matlabbatch{1}.spm.stats.fmri_spec.mask             = {''};
    matlabbatch{1}.spm.stats.fmri_spec.cvi              = 'AR(1)';
    matlabbatch{1}.spm.stats.fmri_spec.dir              = {[dir_root,'bids_data/derivatives/',subj_handle,'/spm_denoised']};
    matlabbatch{1}.spm.stats.fmri_spec.fact             = struct('name', {}, 'levels', {});
    matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi       = {''};
    matlabbatch{1}.spm.stats.fmri_spec.sess.regress     = struct('name', {}, 'val', {});
    matlabbatch{1}.spm.stats.fmri_spec.sess.hpf         = 128;
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans       = all_scans;
    matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg   = {''};   
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
        
    % specify model
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
    % estimate GLM
    matlabbatch{1}.spm.stats.fmri_est.write_residuals           = 0;
    matlabbatch{1}.spm.stats.fmri_est.method.Classical          = 1;
    matlabbatch{1}.spm.stats.fmri_est.spmmat                    = {[dir_root,'bids_data/derivatives/',subj_handle,'/spm_denoised/SPM.mat']};
      
    % contrast betas
    matlabbatch{2}.spm.stats.con.delete                         = 0;
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.name           = 'contrastModality_atEncoding';
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.convec         = [1 1 -1 -1 0 0 0 0];
    matlabbatch{2}.spm.stats.con.consess{1}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.name           = 'contrastModality_atRetrieval';
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.convec         = [0 0 0 0 1 1 -1 -1];
    matlabbatch{2}.spm.stats.con.consess{2}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{3}.tcon.name           = 'retrievalSuccess_visual';
    matlabbatch{2}.spm.stats.con.consess{3}.tcon.convec         = [0 0 0 0 1 -1 0 0];
    matlabbatch{2}.spm.stats.con.consess{3}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.consess{4}.tcon.name           = 'retrievalSuccess_auditory';
    matlabbatch{2}.spm.stats.con.consess{4}.tcon.convec         = [0 0 0 0 0 0 1 -1];
    matlabbatch{2}.spm.stats.con.consess{4}.tcon.sessrep        = 'none';
    matlabbatch{2}.spm.stats.con.spmmat(1)                      = {[dir_root,'bids_data/derivatives/',subj_handle,'/spm_denoised/SPM.mat']};    
    
    % run job
    spm_jobman('run',matlabbatch)
    clear matlabbatch
        
end

%% Run Second-Level Stats
% define each first-level contrast name
contrast_labels = {'contrastModality_atEncoding' 'contrastModality_atRetrieval' 'retrievalSuccess_visual' 'retrievalSuccess_auditory'};

% cycle through each first-level contrast
for i = 1 : numel(contrast_labels)
    
    % define directory for second-level analysis
    stat_dir = [dir_root,'bids_data/derivatives/group/spm_denoised/',contrast_labels{i},'/'];
    mkdir(stat_dir)
    
    % predefine cell to hold contrast filenames
    contrast_files = cell(n_subj,1);
    
    % cycle through each subkect
    for subj = 1 : n_subj
        
        % add filename of contrast image to group cell
        subj_handle = sprintf('sub-%02.0f',subj);
        contrast_files{subj,1} = [dir_root,'bids_data/derivatives/',subj_handle,'/spm_denoised/con_',sprintf('%04.0f',i),'.nii'];
        
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