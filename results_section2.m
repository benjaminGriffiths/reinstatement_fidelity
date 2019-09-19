%% Power Analysis
% clear workspace
clearvars
close all
clc

% define root directory
if ispc;    dir_bids = 'Y:/projects/reinstatement_fidelity/bids_data/';
            dir_tool = 'Y:/projects/general/';
            dir_repos = 'E:/bjg335/projects/'; % repository directory
else;       dir_bids = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/bids_data/';
            dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';
            dir_repos = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/git_clone/'; % repository directory
end

% add subfunctions
addpath([dir_repos,'reinstatement_fidelity/subfunctions'])
addpath(genpath([dir_repos,'fooof_mat/']))

% define number of subjects
n_subj = 21;

%% Get Power and Slope for Each Participant
% cycle through each subject
tic
for subj = 1 : n_subj
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];    
    
    % load in raw data
    fprintf('\nloading sub-%02.0f data...\n',subj);
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'])  
    
    % run fooof analysis
    [freq,spec,h] = get_fooof(source,'visual');
    
    % save data
    save([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-fooof.mat'],'freq','spec')  
    saveas(h,sprintf('%sderivatives/group/figures/sub-%02.0f_plot-fooof.jpg',dir_bids,subj),'jpg')  
    
    % update command line
    fprintf('\nsub-%02.0f complete...\ntotal time elapsed: %1.1f hours...\nestimated time remaining: %1.1f hours...\n',subj,round(toc/3600,1),round(((toc/subj)*n_subj)/3600,1))
    
    % tidy up
    close all
    clear dir_subj subj freq spec h source
end
    
% save data
mkdir([dir_bids,'derivatives/group/eeg/'])
save([dir_bids,'derivatives/group/eeg/group_task-all_eeg-group_fooof.mat'],'group_freq','-v7.3');

%% Run Per-Participant Regression
% predefine group structure
group_betas = cell(n_subj,4);

% cycle through each subject
for subj = 1 : 14% n_subj
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];    
    
    % load in fooof data
    fprintf('\nloading sub-%02.0f data...\n',subj);
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-fooof.mat']) 
    
    % run fooof analysis
    [group_betas(subj,:),h] = get_betas(freq);
    
    % save figure    
    saveas(h,sprintf('%sderivatives/group/figures/sub-%02.0f_plot-regression.jpg',dir_bids,subj),'jpg')  
    close all
    clear h spec freq subj dir_subj
end
    
%% Get Group Average
% load group data if needed
if ~exist('group_irasa','var')
    load([dir_bids,'derivatives/group/eeg/group_task-all_eeg-group_irasa.mat'])
end

% get number of channels
n_chan = numel(group_irasa{1,1,1}.label);

% create cell for peak frequencies
peak_idx = nan(n_subj,n_chan); tic
                 
% cycle through subjects
for subj = 1 : n_subj
        
    % cycle through channels
    for chan = 1 : n_chan

        % get average oscillatory spectrum
        osci_spec = cat(1,group_irasa{subj,1,1}.prepow(:,chan,:),group_irasa{subj,1,1}.powspctrm(:,chan,:),group_irasa{subj,2,1}.prepow(:,chan,:),group_irasa{subj,2,1}.powspctrm(:,chan,:));
        osci_spec = squeeze(mean(osci_spec,1));

        % find peaks of spectrum
        [vals,idx] = findpeaks(osci_spec);

        % if peaks are present
        if ~isempty(idx)
            [~,midx] = max(vals);
            peak_idx(subj,chan) = idx(midx);
        end
    end

    % update user
    fprintf('peaks estimated for sub-%02.0f... time elapsed: %1.1f minutes\n',subj,toc/60);
end
     
% predefine cells to house concatenated group data
grand_irasa = repmat({struct('cfg',[],'freq',8,'time',1,'label',{group_irasa{1,1,1}.label},...
                     'dimord','subj_chan_freq_time','powspctrm',nan(n_subj,n_chan),...
                     'powslope',nan(n_subj,n_chan),'powoffset',nan(n_subj,n_chan))},[5 1]); tic          
        
% cycle through subjects
for subj = 1 : n_subj
        
    % cycle through channels
    for chan = 1 : n_chan

        % define linear name
        ln = {'trl_at_encoding','trl_at_retrieval'};
        
        % cycle through encoding and retrieval
        for i = 1 : 2

            % get number of trials        
            n_trl  = numel(group_irasa{subj,i,1}.trialinfo);
        
            % create pre/post model
            Y = [ones(n_trl,1);zeros(n_trl,1)]+1;
            
            % create regressor table
            X        = nan(n_trl*2,3);
            X(:,1)   = cat(1,group_irasa{subj,i,1}.powslope(:,chan),group_irasa{subj,i,1}.preslope(:,chan));
            X(:,2)   = cat(1,group_irasa{subj,i,1}.powoffset(:,chan),group_irasa{subj,i,1}.preoffset(:,chan));
            
            % get power regressor if peak found
            idx = peak_idx(subj,chan);
            if ~isnan(idx)
                X(:,3)   = mean(cat(1,group_irasa{subj,i,1}.powspctrm(:,chan,idx-1:idx+1),group_irasa{subj,i,1}.prepow(:,chan,idx-1:idx+1)),3);
            end
            
            % add linear regressor
            for trl = 1 : n_trl
                X(trl,4) = group_irasa{subj,i,1}.trialinfo{trl}.(ln{i});
            end
            
            % handle NaNs
            nidx = any(isnan(X),2);
            X(nidx,:) = [];
            Y(nidx,:) = [];
            
            % run multiregression (if peak found)
            if ~isempty(X)
                B = fitlm(X,Y);

                % store data
                grand_irasa{i,1}.powslope(subj,chan)    = B.Coefficients.tStat(2);
                grand_irasa{i,1}.powoffset(subj,chan)   = B.Coefficients.tStat(3);
                grand_irasa{i,1}.powspctrm(subj,chan)   = B.Coefficients.tStat(4);
            end
        end           
    end    
    
    % update user
    fprintf('sub-%02.0f complete... time elapsed: %1.1f minutes\n',subj,round(toc/60,1));
    clear j k n_trl tmp_pd tmp_pow tmp_pp mid pid vals idx valp idxp erd    
end

% cycle through subjects
for subj = 1 : n_subj

    % get data of interest
    data = group_irasa{subj,2,1};
    
    % get number of trials        
    n_trl  = numel(data.trialinfo);
    
    % cycle through channels
    for chan = 1 : n_chan
        
        % create regressor table
        X        = nan(n_trl,7);
        X(:,1)   = data.powslope(:,chan);
        X(:,2)   = data.powoffset(:,chan);
        X(:,4)   = data.preslope(:,chan);
        X(:,5)   = data.preoffset(:,chan);
                        
        % get power regressor if peak found
        idx = peak_idx(subj,chan);
        if ~isnan(idx)
            X(:,3)   = mean(data.powspctrm(:,chan,idx-1:idx+1),3);
            X(:,6)   = mean(data.prepow(:,chan,idx-1:idx+1),3);
        end

        % create memory model
        Y = nan(n_trl,1);
        for i = 1 : n_trl
            Y(i) = data.trialinfo{i}.recalled;
            X(i,7) = data.trialinfo{i}.trl_at_retrieval;
        end    
        
        % handle NaNs
        idx = any(isnan(X),2);
        X(idx,:) = [];        
        Y(idx,:) = [];

        % run multiregression (if peak found)
        if ~isempty(X)
            
            % run multiple regression
            B = fitlm(X,Y);

            % store poststim data
            grand_irasa{3,1}.powslope(subj,chan)    = B.Coefficients.tStat(2);
            grand_irasa{3,1}.powoffset(subj,chan)   = B.Coefficients.tStat(3);
            grand_irasa{3,1}.powspctrm(subj,chan)   = B.Coefficients.tStat(4);
            
            % store prestim data
            grand_irasa{4,1}.powslope(subj,chan)    = B.Coefficients.tStat(5);
            grand_irasa{4,1}.powoffset(subj,chan)   = B.Coefficients.tStat(6);
            grand_irasa{4,1}.powspctrm(subj,chan)   = B.Coefficients.tStat(7);
            
            % create difference model
            Xd = zeros(size(X,1),4);
            Xd(:,1) = X(:,1) - X(:,4);
            Xd(:,2) = X(:,2) - X(:,5);
            Xd(:,3) = X(:,3) - X(:,6);
            Xd(:,4) = X(:,7);
            
            % run multiple regression
            B = fitlm(X,Y);

            % store difference data
            grand_irasa{5,1}.powslope(subj,chan)    = B.Coefficients.tStat(2);
            grand_irasa{5,1}.powoffset(subj,chan)   = B.Coefficients.tStat(3);
            grand_irasa{5,1}.powspctrm(subj,chan)   = B.Coefficients.tStat(4);           
        end
    end

    % update user
    fprintf('sub-%02.0f complete... time elapsed: %1.1f minutes\n',subj,toc/60);
end
    
% save data
save([dir_bids,'derivatives/group/eeg/group_task-all_eeg-irasa.mat'],'grand_irasa');

%% Change Datatype to Source
% load sourcemodel
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat']);

% load AAL atlas
mri = ft_read_mri([dir_bids,'sourcedata/masks/whole_brain.nii']);
 
% interpolate clusters with grid
cfg             = [];
cfg.parameter	= 'anatomy';
roi             = ft_sourceinterpolate(cfg,mri,sourcemodel);

% define additional roi parameters
roi.inside      = ~isnan(roi.anatomy) & roi.anatomy > 0;
roi.anatomy     = double(~isnan(roi.anatomy) & roi.anatomy > 0);
roi.pos         = sourcemodel.pos;

% cycle through each data type
for i = 1 : numel(grand_irasa)
    
    % add in source model details
    grand_irasa{i,1}.powdimord   = 'rpt_pos';
    grand_irasa{i,1}.dim         = roi.dim;
    grand_irasa{i,1}.inside      = roi.inside(:);
    grand_irasa{i,1}.pos         = roi.pos;
    grand_irasa{i,1}.pow         = zeros(n_subj,numel(grand_irasa{i,1}.inside));
    grand_irasa{i,1}.cfg         = [];
    
    % remove freq details
    grand_irasa{i,1} = rmfield(grand_irasa{i,1},{'label','freq','time','dimord'});
end

% create source data strucutre
grand_freq = repmat(grand_irasa,[1 3]);
powlabel = {'powspctrm','powslope','powoffset'};

% cycle through each condition and metric
for i = 1 : size(grand_freq,1)
    for j = 1 : size(grand_freq,2)

        % add in "outside" virtual electrods
        tmp_pow = zeros(size(grand_freq{i,j}.pow,1),size(grand_freq{i,j}.inside,1));
        tmp_pow(:,grand_freq{i,j}.inside) = grand_freq{i,j}.(powlabel{j});
        grand_freq{i,j}.pow = tmp_pow;
        
        % remove irrelevant fields
        grand_freq{i,j} = rmfield(grand_freq{i,j},powlabel);
    end
end

%% Run Statistics
% define cluster names
plot_names = {'enc_pow','ret_pow','rse_pow','pre_pow','diff_pow',...
    'enc_slope','ret_slope','rse_slope','pre_slope','diff_slope',...
    'enc_offset','ret_offset','rse_offset','pre_offset','diff_offset'};

% set seed
rng(1) 

% test oscillatory power
cfg             = [];
cfg.tail        = [zeros(3,1)-1; 1; -1; zeros(10,1)];
cfg.parameter   = 'pow';
[stat,tbl]      = run_oneSampleT(cfg, grand_freq(:));
tbl             = addvars(tbl,plot_names','Before','t');

% save data
save([dir_bids,'derivatives/group/eeg/group_task-all_pow-stat.mat'],'stat','tbl');

%% Extract Raw Power of Cluster 
% prepare table for stat values
tbl = array2table(zeros(n_subj,numel(plot_names)),'VariableNames',plot_names);

% cycle through conditions
for i = 1 : 15
     
    % get indices of clusters
    if i > 5 || i == 4
        clus_idx = stat{i}.posclusterslabelmat==1;
    else
        clus_idx = stat{i}.negclusterslabelmat==1;
    end

    % create table
    tbl.(plot_names{i})(:,1) = nanmean(grand_freq{i}.pow(:,clus_idx),2);
end

% write table
writetable(tbl,[dir_bids,'derivatives/group/eeg/group_task-all_eeg-cluster.csv'],'Delimiter',',')
copyfile([dir_bids,'derivatives/group/eeg/group_task-all_eeg-cluster.csv'],[dir_repos,'data/fig2_data/group_task-rf_eeg-cluster.csv'])

%% Create Masked Surface Plots
% load template MRI
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat']);

% cycle through conditions
for i = 1 : numel(plot_names)
        
    % get indices of clusters
    if i > 5 || i == 4
        clus_idx = stat{i}.posclusterslabelmat==1;
    else
        clus_idx = stat{i}.negclusterslabelmat==1;
    end
    
    % create source data structure
    source                  = [];
    source.inside           = stat{i}.inside;
    source.dim              = stat{i}.dim;
    source.pos              = stat{i}.pos*10;
    source.unit             = 'mm';
    
    % define powspctrm of cluster
    source.pow              = nan(size(stat{i}.pos,1),1);     
    source.pow(clus_idx)	= stat{i}.stat(clus_idx); 
    
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
    cfg.filename      = [dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-maskedmap.nii'];  % enter the desired file name
    cfg.filetype      = 'nifti';
    cfg.coordsys      = 'spm';
    cfg.vmpversion    = 2;
    cfg.datatype      = 'float';
    ft_volumewrite(cfg,source);      % be sure to use your interpolated source data
    
    % reslice to 1mm isometric to match template MRI
    reslice_nii([dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-maskedmap.nii'],[dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-maskedmap.nii'],[1 1 1]);
    copyfile([dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-maskedmap.nii'],[dir_repos,'data/fig2_data/group_task-',plot_names{i},'_eeg-maskedmap.nii'])    
end

%% Create Unasked Surface Plots
% load template MRI
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat']);

% cycle through conditions
for i = 1 : numel(plot_names)
      
    % create source data structure
    source                  = [];
    source.inside           = stat{i}.inside;
    source.dim              = stat{i}.dim;
    source.pos              = stat{i}.pos*10;
    source.unit             = 'mm';
    
    % define powspctrm of cluster
    source.pow              = nan(size(stat{i}.pos,1),1);     
    source.pow(source.inside) = stat{i}.stat(source.inside); 
    
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
    cfg.filename      = [dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-unmaskedmap.nii'];  % enter the desired file name
    cfg.filetype      = 'nifti';
    cfg.coordsys      = 'spm';
    cfg.vmpversion    = 2;
    cfg.datatype      = 'float';
    ft_volumewrite(cfg,source);      % be sure to use your interpolated source data
    
    % reslice to 1mm isometric to match template MRI
    reslice_nii([dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-unmaskedmap.nii'],[dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-unmaskedmap.nii'],[1 1 1]);
    copyfile([dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-unmaskedmap.nii'],[dir_repos,'data/fig2_data/group_task-',plot_names{i},'_eeg-unmaskedmap.nii'])    
end
