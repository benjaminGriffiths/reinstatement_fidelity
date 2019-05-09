%% Behavioural
% initialisation
clearvars
clc

% define root directory
if ispc     
    dir_root = 'Y:/projects/reinstatement_fidelity/';
    dir_bids = [dir_root,'bids_data/'];
    dir_tool = 'Y:/projects/general/';
    dir_repos = 'E:/bjg335/projects/reinstatement_fidelity/'; % repository directory
    
    % add subfunctions
    addpath([dir_repos,'subfunctions'])
    
else       
    dir_root = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/';
    dir_bids = [dir_root,'bids_data/'];
    dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';            
    dir_repos = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/scripts/'; % repository directory
    
    % add subfunctions
    addpath([dir_repos,'subfunctions'])
end

% cycle through each subject
for subj = 1 : 21
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    
    % create mat to hold memory data
    mem = [];
        
    % cycle through each encoding run
    for run = [1 3 5 7]

        % load event table
        tbl = readtable([dir_root,'bids_data/',subj_handle,'/func/',subj_handle,'_task-rf_run-',num2str(run),'_events.tsv'],'FileType','text','Delimiter','\t');
        
        % cut table to stimulus onset
        tbl = tbl(ismember(tbl.trial_type,'Stimulus Onset'),:);
        
        % extract memory
        mem(end+1:end+size(tbl,1)) = tbl.recalled;
    end
    
    % get mean
    avg_mem(subj) = nanmean(mem);
end

