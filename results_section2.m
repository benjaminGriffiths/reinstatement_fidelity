%% Power Analysis
% clear workspace
clearvars
close all
clc

% start parallel pool
p = gcp('nocreate');
if isempty(p)
    parpool(4);
end

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
% predefine group freq structure
group_freq = cell(n_subj,6); tic

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];    
    
    % load in raw data
    fprintf('\nloading sub-%02.0f data...\n',subj);
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'])  
    
    % run fooof analysis
    [freq,sigout,h] = get_condition_fooof(source);
    
    % save outputs
    saveas(h,sprintf('%sderivatives/group/figures/fooof/sub-%02.0f_task-vis_eeg-fooof.jpg',dir_bids,subj),'jpg')
    save(sprintf('%sderivatives/sub-%02.0f/eeg/sub-%02.0f_task-rf_eeg-fooof.mat',dir_bids,subj,subj),'freq')  
    
    % add data to group
    group_freq(subj,:) = freq;
    
    % extract signal
    signal(subj,:,1) = mean(sigout{1,1}(:,3,:),3);
    signal(subj,:,2) = mean(sigout{1,1}(:,4,:),3);
    signal(subj,:,3) = mean(sigout{1,1}(:,5,1),3);
    
    % update command line
    te = round(toc/60,1);
    tr = round((te/subj)*(n_subj-subj),1);
    fprintf('\nsub-%02.0f complete...\ntotal time elapsed: %1.1f minutes...\nestimated time remaining: %1.1f minutes...\n',subj,te,tr)
    
    % tidy up
    close all
    clear dir_subj subj freq spec h source te tr
end

% plot outcome
figure; hold on
f = squeeze(signal(1,:,3));
shadedErrorBar(f,mean(signal(:,:,1)),sem(signal(:,:,1)))
shadedErrorBar(f,mean(signal(:,:,2)),sem(signal(:,:,2)))

%% Get Group Average
% create grand data cell
grand_freq = cell(size(group_freq,2),1);

% cycle through each condition
for i = 1 : size(group_freq,2)
    
    % get grand average
    cfg                 = [];
    cfg.keepindividual  = 'yes';
    cfg.parameter       = {'powspctrm','slope'};
    grand_freq{i,1}     = ft_freqgrandaverage(cfg,group_freq{:,i});  
end

% tidy workspace
clear i group_betas cfg

%% Change Datatype to Source
% load ROI data
load([dir_bids,'derivatives/group/eeg/roi.mat'])

% cycle through each data type
for i = 1 : numel(grand_freq)
    
    % add in source model details
    grand_freq{i,1}.powdimord   = 'rpt_pos';
    grand_freq{i,1}.dim         = roi.dim;
    grand_freq{i,1}.inside      = roi.inside(:);
    grand_freq{i,1}.pos         = roi.pos;
    grand_freq{i,1}.pow         = zeros(n_subj,numel(grand_freq{i,1}.inside));
    grand_freq{i,1}.slp         = zeros(n_subj,numel(grand_freq{i,1}.inside));
    grand_freq{i,1}.cfg         = [];
    
    % add in "outside" virtual electrodes to pow
    tmp_pow = zeros(size(grand_freq{i,1}.pow,1),size(grand_freq{i,1}.inside,1));
    tmp_pow(:,grand_freq{i,1}.inside) = grand_freq{i,1}.powspctrm;
    grand_freq{i,1}.pow = tmp_pow;
        
    % add in "outside" virtual electrodes to slope
    tmp_pow = zeros(size(grand_freq{i,1}.pow,1),size(grand_freq{i,1}.inside,1));
    tmp_pow(:,grand_freq{i,1}.inside) = grand_freq{i,1}.slope;
    grand_freq{i,1}.slp = tmp_pow;
    
    % remove freq details
    grand_freq{i,1} = rmfield(grand_freq{i,1},{'powspctrm','slope','label','freq','time','dimord','cfg'});
    clear tmp_pow i
end

% save data
save([dir_bids,'derivatives/group/eeg/group_task-all_pow-fooof.mat'],'grand_freq');

%% Run Statistics
% set seed
rng(1) 

% test oscillatory power
cfg             = [];
cfg.tail        = zeros(6,1)-1;
cfg.parameter   = 'pow';
[s_pow,t_pow]   = run_oneSampleT(cfg, grand_freq(:));

% test slope
cfg             = [];
cfg.tail        = zeros(6,1);
cfg.parameter   = 'slp';
[s_slp,t_slp]   = run_oneSampleT(cfg, grand_freq(:));

% combine tables
stat = cat(1,s_pow,s_slp);
tbl = cat(1,t_pow,t_slp);

% add labels
modality = repmat({'visual','audio'},[6 1]);
measure  = repmat({'power','slope'},[3 1]);
task     = repmat({'encoding','retrieval','RSE'},[1 4]);
tbl = addvars(tbl,modality(:),'Before','t');
tbl = addvars(tbl,repmat(measure(:),[2 1]),'Before','t');
tbl = addvars(tbl,task(:),'Before','t');
disp(tbl)

% save data
save([dir_bids,'derivatives/group/eeg/group_task-all_pow-stat.mat'],'stat','tbl');

%% Extract Raw Power of Cluster 
% prepare table for stat values
plot_names={'enc','ret','rse'};
tbl = array2table(zeros(n_subj,numel(plot_names)),'VariableNames',plot_names);
% cycle through conditions
for i = 1 : 3%numel(plot_names)
     
    % get indices of power clusters
    clus_idx = stat{i}.negclusterslabelmat==1;

    % create table
    tbl.(plot_names{i})(:,1) = nanmean(grand_freq{i}.pow(:,clus_idx),2);
    
    % get indices of slope clusters
    if i ~= 4
        clus_idx = stat{i}.negclusterslabelmat==1;
    else
        clus_idx = stat{i}.posclusterslabelmat==1;
    end

    % create table
    tbl.(plot_names{i})(:,1) = nanmean(grand_freq{i}.slp(:,clus_idx),2);
end

% get rough plot
figure; hold on
for i = 1 : 3
    plot(i,tbl.(plot_names{i}),'ko')
end
set(gca,'xtick',1:3,'xticklabel',plot_names)
xlim([0 4])

% write table
writetable(tbl,[dir_bids,'derivatives/group/eeg/group_task-all_eeg-cluster.csv'],'Delimiter',',')
copyfile([dir_bids,'derivatives/group/eeg/group_task-all_eeg-cluster.csv'],[dir_repos,'reinstatement_fidelity/data/fig2_data/group_task-rf_eeg-cluster.csv'])

%% Create Masked Surface Plots
% load template MRI
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat']);

% cycle through conditions
for i = 1 :3% numel(plot_names)
        
    % get indices of clusters
    if i ~= 6
        clus_idx = stat{i}.negclusterslabelmat==1;
    else
        clus_idx = stat{i}.posclusterslabelmat==1;
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
    copyfile([dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-maskedmap.nii'],[dir_repos,'reinstatement_fidelity/data/fig2_data/group_task-',plot_names{i},'_eeg-maskedmap.nii'])    
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
