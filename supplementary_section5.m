%% Phase Analysis
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
mri = ft_read_mri([dir_bids,'sourcedata/masks/whole_brain.nii']);
 
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
%group_freq = cell(n_subj,2);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];
    
    % load in raw data
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'])  
    
    % get time-frequency representaiton of data for specified conditions
    group_freq{subj,1} = get_phasepowercorr(source,'encoding','visual');
    %group_freq{subj,2} = get_phasepowercorr(source,'retrieval','visual');
    
    % update command line
    fprintf('\nSubject %02.0f of %02.0f complete...\n\n',subj,n_subj)
end
    
% predefine cells to house concatenated group data
grand_freq = cell(size(group_freq,2),1);

% cycle through conditions
for i = 1 : size(group_freq,2)
    
    % get grand average of subjects
    cfg                 = [];
    cfg.keepindividual  = 'yes';
    grand_freq{i,1}     = ft_freqgrandaverage(cfg,group_freq{:,i});
    
    % add in source model details
%     grand_freq{i,1}.powdimord   = 'rpt_pos';
%     grand_freq{i,1}.dim         = roi.dim;
%     grand_freq{i,1}.inside      = roi.inside(:);
%     grand_freq{i,1}.pos         = roi.pos;
%     grand_freq{i,1}.cfg         = [];
%     
%     % add in "outside" virtual electrods
%     tmp_pow = zeros(size(grand_freq{i,1}.powspctrm,1),size(grand_freq{i,1}.inside,1));
%     tmp_pow(:,grand_freq{i,1}.inside) = grand_freq{i,1}.powspctrm;
%     grand_freq{i,1}.pow = tmp_pow;
%     
    % remove freq details
    %grand_freq{i,1} = rmfield(grand_freq{i,1},{'powspctrm','label','freq','time','dimord'});
end
    
% save data
save([dir_bids,'derivatives/group/eeg/group_task-all_eeg-phasepowercorr.mat'],'grand_freq');

%% Run Statistics
% average over roi
load([dir_bids,'derivatives/group/eeg/group_task-all_eeg-phasestat.mat'])

% average roi
grand_freq{1}.powspctrm = nanmean(grand_freq{1}.powspctrm(:,stat{1}.posclusterslabelmat(stat{1}.inside)==1),2);
grand_freq{2}.powspctrm = nanmean(grand_freq{2}.powspctrm(:,stat{1}.posclusterslabelmat(stat{1}.inside)==1),2);
grand_freq{1}.label = {'dummy'};
grand_freq{2}.label = {'dummy'};

% predefine cell for statistics
cfg = [];
cfg.tail            = -1;
%cfg.parameter       = 'pow';
[stat,tbl]          = run_oneSampleT(cfg, grand_freq);

% save data
save([dir_bids,'derivatives/group/eeg/group_task-all_eeg-phasepowerstat.mat'],'stat','tbl');

%% Extract Raw Power of Cluster 
% define cluster names
plot_names = {'eesim','ersim'};

% prepare table for stat values
tbl = array2table(zeros(n_subj,numel(plot_names)),'VariableNames',plot_names);

% cycle through conditions
for i = 1 : numel(stat)
     
    % get indices of clusters
    clus_idx = stat{i}.posclusterslabelmat==1;

    % create table
    tbl.(plot_names{i})(:,1) = nanmean(grand_freq{i}.pow(:,clus_idx),2);
end

% write table
writetable(tbl,[dir_bids,'derivatives/group/eeg/group_task-all_eeg-phasecluster.csv'],'Delimiter',',')
copyfile([dir_bids,'derivatives/group/eeg/group_task-all_eeg-phasecluster.csv'],[dir_repos,'data/sup2_data/group_task-rf_eeg-phasecluster.csv'])

%% Create Surface Plots
% load template MRI
load([dir_tool,'fieldtrip-20170319/template/sourcemodel/standard_sourcemodel3d10mm.mat']);

% define plot names
plot_names = {'eesim','ersim'};

% cycle through conditions
for i = 1 : numel(stat)
        
    % get indices of clusters
    clus_idx = stat{i}.posclusterslabelmat==1;

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
    cfg.filename      = [dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-map.nii'];  % enter the desired file name
    cfg.filetype      = 'nifti';
    cfg.coordsys      = 'spm';
    cfg.vmpversion    = 2;
    cfg.datatype      = 'float';
    ft_volumewrite(cfg,source);      % be sure to use your interpolated source data
    
    % reslice to 1mm isometric to match template MRI
    reslice_nii([dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-map.nii'],[dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-map.nii'],[1 1 1]);
    copyfile([dir_bids,'derivatives/group/eeg/group_task-',plot_names{i},'_eeg-map.nii'],[dir_repos,'data/sup2_data/group_task-',plot_names{i},'_eeg-map.nii'])    
end