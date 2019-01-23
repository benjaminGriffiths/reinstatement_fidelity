%% Power Analysis
% clear workspace
clearvars
close all
clc

% define root directory
if ispc;    dir_bids = 'Y:/projects/reinstatement_fidelity/bids_data/';
            dir_tool = 'Y:/projects/general/';
            dir_repos = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
else;       dir_bids = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/bids_data/';
            dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';
            dir_repos = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/git_clone/reinstatement_fidelity/'; % repository directory
end

% add subfunctions
addpath([dir_repos,'subfunctions'])

% define number of subjects
n_subj = 21;

%% Prepare ROI
% load sourcemodel
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat']);

% load AAL atlas
mri = ft_read_mri([dir_bids,'sourcedata/masks/aal.nii']);
 
% interpolate clusters with grid
cfg             = [];
cfg.parameter	= 'anatomy';
roi             = ft_sourceinterpolate(cfg,mri,sourcemodel);

% define additional roi parameters
roi.inside      = ~isnan(roi.anatomy) & roi.anatomy > 0;
roi.anatomy     = double(~isnan(roi.anatomy) & roi.anatomy > 0);
roi.pos         = sourcemodel.pos;

%% Get Timelocked Representation of Each Condition
% predefine cell for group data
group_freq = cell(n_subj,2);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];
    
    % load in raw data
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'])  
    
    % get time-frequency representaiton of data for specified conditions
    group_freq{subj,1} = get_basediff_timefrequency(source,'encoding','visual');
    group_freq{subj,2} = get_memdiff_timefrequency(source,'retrieval','visual');
end
    
% predefine cells to house concatenated group data
grand_freq = cell(size(group_freq,2),1);

% cycle through conditions
for i = 1 : size(group_freq,2)
    
    % get grand average of subjects
    cfg                 = [];
    cfg.keepindividual  = 'yes';
    grand_freq{i,1}     = ft_freqgrandaverage(cfg,group_freq{:,i});
    
    % collapse over time/frequency
    grand_freq{i,1}.pow = mean(mean(grand_freq{i,1}.powspctrm,4),3);
    
    % add in source model details
    grand_freq{i,1}.powdimord   = 'rpt_pos';
    grand_freq{i,1}.dim         = roi.dim;
    grand_freq{i,1}.inside      = roi.inside(:);
    grand_freq{i,1}.pos         = roi.pos;
    grand_freq{i,1}.cfg         = [];
    
    % add in "outside" virtual electrods
    tmp_pow = zeros(size(grand_freq{i,1}.pow,1),size(grand_freq{i,1}.inside,1));
    tmp_pow(:,grand_freq{i,1}.inside) = grand_freq{i,1}.pow;
    grand_freq{i,1}.pow = tmp_pow;
    
    % remove freq details
    grand_freq{i,1} = rmfield(grand_freq{i,1},{'powspctrm','label','freq','time','dimord'});
end
    
% save data
save([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-source.mat'],'grand_freq');

%% Run Statistics
% predefine cell for statistics
cfg.tail        = -1;
cfg.parameter   = 'pow';
[stat,tbl]      = run_oneSampleT(cfg, grand_freq);

% save data
save([dir_bids,'derivatives/group/eeg/group_task-all_eeg-stat.mat'],'stat','tbl');

%% Extract Raw Power of Cluster 
% prepare table for stat values
tbl = array2table(zeros(n_subj,2),'VariableNames',plot_names);

% cycle through conditions
for i = 1 : numel(stat)
     
    % get indices of clusters
    clus_idx = stat{i}.negclusterslabelmat==1;

    % create table
    tbl.(plot_names{i})(:,1) = nanmean(grand_freq{i}.pow(:,clus_idx),2);
end

% write table
writetable(tbl,[dir_bids,'derivatives/group/eeg/group_task-rf_eeg-cluster.csv'],'Delimiter',',')

%% Get Time/Frequency Series of Statistical Cluster
% predefine cell for group data
group_freq = cell(n_subj,2);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];
    
    % load in raw data
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'])  
        
    % define custom TF resolution
    res.toi = -1:0.05:3;
    res.foi = 3:40;
    
    % restrict data to cluster
    cfg         = [];
    cfg.channel = source.label(stat{1}.negclusterslabelmat(stat{1}.inside)==1);
    sourcetmp   = ft_selectdata(cfg,source);

    % get time-frequency representaiton of data for specified conditions
    group_freq{subj,1} = get_basediff_timefrequency(sourcetmp,'encoding','visual',res);
    
    % restrict data to cluster
    cfg         = [];
    cfg.channel = source.label(stat{2}.negclusterslabelmat(stat{2}.inside)==1);
    sourcetmp   = ft_selectdata(cfg,source);
    
    % get time-frequency representaiton of data for specified conditions
    group_freq{subj,2} = get_memdiff_timefrequency(sourcetmp,'retrieval','visual',res);
end
    
% predefine cells to house concatenated group data
grand_freq = cell(size(group_freq,2),1);

% cycle through conditions
for i = 1 : size(group_freq,2)
    
    % get grand average of subjects
    cfg                 = [];
    cfg.keepindividual  = 'yes';
    grand_freq{i,1}     = ft_freqgrandaverage(cfg,group_freq{:,i});
end

% extract encoding/retrieval time/frequency series
enc_freq = squeeze(mean(grand_freq{1}.powspctrm,2));
ret_freq = squeeze(mean(mean(grand_freq{2}.powspctrm(:,:,:,grand_freq{2}.time>=0.5 & grand_freq{2}.time<=1.5),4),2));
ret_time = squeeze(nanmean(nanmean(grand_freq{2}.powspctrm(:,:,grand_freq{2}.freq>=8 & grand_freq{2}.freq<=30,:),3),2));

% save details
csvwrite([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-encfreqspc.csv'],enc_freq)
csvwrite([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-retfreqspc.csv'],ret_freq)
csvwrite([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-rettimespc.csv'],ret_time)

%% Create Surface Plots
% load template MRI
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat']);

% define plot names
plot_names = {'encoding','retrieval'};

% cycle through conditions
for i = 1 : numel(stat)
        
    % get indices of clusters
    clus_idx = stat{i}.negclusterslabelmat==1;

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
    cfg.filename      = [dir_bids,'derivatives/group/eeg/group_task-rf_eeg-',plot_names{i},'map.nii'];  % enter the desired file name
    cfg.filetype      = 'nifti';
    cfg.coordsys      = 'spm';
    cfg.vmpversion    = 2;
    cfg.datatype      = 'float';
    ft_volumewrite(cfg,source);      % be sure to use your interpolated source data
    
    % reslice to 1mm isometric to match template MRI
    reslice_nii([dir_bids,'derivatives/group/eeg/group_task-rf_eeg-',plot_names{i},'map.nii'],[dir_bids,'derivatives/group/eeg/group_task-rf_eeg-',plot_names{i},'map.nii'],[1 1 1]);    
end

