%% RSA
% initialisation
clearvars
clc

% initalise SPM
spm('Defaults','fMRI');
spm_jobman('initcfg');

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
end

% add RSA toolbox to path
addpath(genpath([dir_tool,'rsatoolbox-develop']))

%% Define Key Parameters
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

%% Prepare GLM and Data for Searchlight Analysis
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/derivatives/',subj_handle,'/'];
    
    % load pattern data
    load([dir_subj,'/rsa/',subj_handle,'_task-rf_rsa-maskedDemeanedVolume.mat'])
    
    % load SPM.mat
    load([dir_subj,'/rsa/SPM.mat'])
    
    % create table to record stimulus detail
    stim_details  = array2table(zeros(n_trials*2,2),'VariableNames', {'encoding','modality'});
    stim_count = 1;
    
    % create table to record scan details
    scan_details  = array2table(zeros(n_volumes*8,2),'VariableNames', {'encoding','modality'});
    scan_count = 1;  
    
    % cycle through each run
    for run = 1 : n_runs

        % load event table
        tbl = readtable([dir_root,'bids_data/',subj_handle,'/func/',subj_handle,'_task-rf_run-',num2str(run),'_events.tsv'],'FileType','text','Delimiter','\t');
        
        % check block types
        block_encoding  = double(any(ismember(tbl.operation,'encoding')));
        block_visual    = double(any(ismember(tbl.modality,'Visual')));
        
        % cycle through every event
        for e = 1 : size(tbl,1)
            
            % if an event
            if strcmpi(tbl.trial_type(e),'Stimulus Onset')
            
                % add key values
                stim_details.encoding(stim_count) = strcmpi(tbl.operation(e),'encoding');
                stim_details.modality(stim_count) = strcmpi(tbl.modality(e),'Visual');
                
                % get stimulus values
                switch tbl.stimulus{e}
                    case 'WATERMILL';   stim_details.stimulus(stim_count,1) = 1;
                    case 'UNDERWATER';  stim_details.stimulus(stim_count,1) = 2;
                    case 'BIKE';        stim_details.stimulus(stim_count,1) = 3;
                    case 'FARM';        stim_details.stimulus(stim_count,1) = 4;
                end
                
                % change value if retrieval
                if stim_details.encoding(stim_count) ~= 1 && stim_details.modality(stim_count) == 1
                    stim_details.stimulus(stim_count,1) = stim_details.stimulus(stim_count,1) + 4;
                end
                
                % update counter
                stim_count = stim_count + 1; 
                
            % else if a volume
            elseif strncmpi(tbl.trial_type(e),'Volume',6)
                
                % add key values
                scan_details.encoding(scan_count) = block_encoding;
                scan_details.modality(scan_count) = block_visual;

                % update counter
                scan_count = scan_count + 1;                 
            end          
        end
    end
    
    % kick out first three/last five scans of each run
    scan_details([1:3 251:258 506:513 761:768 1016:1023 1271:1278 1526:1533 1781:1788 2036:2040],:) = [];
    
    % clean up
    clear run tbl e stim_count
    
    % get design matrix (X) and split into two groups (Xa and Xb)
    X.raw = SPM.xX.X;
        
    % remove scans/regressors that are not visual (to
    % computationally demanding to anything more than this)
    X.raw = X.raw(scan_details.modality==1,:);
    
    % split GLM into two groups (train [encoding] and test [retrieval] data)
    X.a = X.raw(1:size(X.raw,1)/2,:);
    X.b = X.raw((size(X.raw,1)/2)+1:end,:);
    
    % --- prepare patterns --- %
    % get activation matrix (Y) and 
    Y.raw = patterns(scan_details.modality==1,:);
    
    % split into two groups (Ya and Yb)
    Y.a = Y.raw(1:size(Y.raw,1)/2,:);
    Y.b = Y.raw(size(Y.raw,1)/2+1:end,:);
        
    % --- remove singular dimensions --- %
    % find zero-value rows
    X.a_badRow = all(X.a(:,1:8)==0,2);
    X.b_badRow = all(X.b(:,1:8)==0,2);
    
    % remove zero-value rows
    X.a(X.a_badRow,:) = [];
    X.b(X.a_badRow,:) = [];
    
    % remove zero-value rows
    Y.a(X.a_badRow,:) = [];
    Y.b(X.b_badRow,:) = [];
        
    % cycle through A/B avg/trl combinations
    fn = {'a','b'};
    for i = 1 : 2
        for j = 1 : size(X.(fn{i}),2)
            
            % if column is singular
            if numel(unique(X.(fn{i})(:,j))) == 1
                X.([fn{i},'_badCol'])(j,1) = true;
            else
                X.([fn{i},'_badCol'])(j,1) = false;
            end
        end
    end
    
    % remove single-valued columns
    X.a(:,X.a_badCol) = [];
    X.b(:,X.b_badCol) = [];
    
    % save data
    save([dir_subj,'/rsa/',subj_handle,'_task-rf_rsa-FormattedVolume.mat'],'X','Y')
    
    clear X Y i j fn patterns subj_handle dir_subj
    
    fprintf('Subject %02.0f of %02.0f prepared...\n',subj,n_subj)
end

%% Define Searchlight Characteristics
% get the radius of searchlight in voxels
voxRadInSearchLight = scan_search ./ scan_vox;

% define distance from searchlight centre to perimeter in voxels
dist2Perimeter = ceil(voxRadInSearchLight);

% create boolean searchlight sphere
[x,y,z] = meshgrid(-dist2Perimeter(1):dist2Perimeter(1),-dist2Perimeter(2):dist2Perimeter(2),-dist2Perimeter(3):dist2Perimeter(3));
sphere  = ((x*scan_vox(1)).^2+(y*scan_vox(2)).^2+(z*scan_vox(3)).^2)<=(scan_search^2);

% clean up
clear x y z voxRadInSearchLight

%% Run Searchlight LDt Analysis
tic

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % load data
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-FormattedVolume.mat'])
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-mask.mat'])
    
    % predefine model rdm
    model_rdm{1}        	= nan(8,8);

    % set items that belong to the different category to 1
    model_rdm{1}(5:8,1:4)	= 1;

    % set items that belong to a same category to -1
    model_rdm{1}(5,1)       = -1;
    model_rdm{1}(6,2)       = -1;
    model_rdm{1}(7,3)       = -1;
    model_rdm{1}(8,4)       = -1;
    
    % set diagonal to NaN
    model_rdm{1}(logical(eye(size(model_rdm{1},1)))) = 0; 

    % define empty brain map
    rdmBrain = zeros(scan_fov);
    
    % predefine voxels to be analysed
    sl_vox = cell(sum(mask_idx),1);
    
    % predefine vector of searchlights to be used
    goodSL = true(size(sl_vox));
    
    % get index of all voxels within mask
    M = find(mask_idx);
    
    % update command line
    fprintf('Defining searchlights...\n')
    
    % run LDt analysis for every voxel
    for vox = 1 : numel(M)
        
        % get subscript co-ordinates of searchlight centre
        [x,y,z] = ind2sub(scan_fov,M(vox));
        
        % define box which houses spherical searchlight
        tmpX = x-dist2Perimeter(1):x+dist2Perimeter(1);
        tmpY = y-dist2Perimeter(2):y+dist2Perimeter(2);
        tmpZ = z-dist2Perimeter(3):z+dist2Perimeter(3);
        
        % remove co-ordinates less than or equal to zero, and greater than FOV
        tmpX(tmpX<=0) = []; tmpX(tmpX>scan_fov(1)) = [];
        tmpY(tmpY<=0) = []; tmpY(tmpY>scan_fov(2)) = [];
        tmpZ(tmpZ<=0) = []; tmpZ(tmpZ>scan_fov(3)) = [];
        
        % get indices of every voxel in sphere
        voxIdx = [];
        for xi = 1 : numel(tmpX)
            for yi = 1 : numel(tmpY)
                for zi = 1 : numel(tmpZ)
                    if sphere(xi,yi,zi)
                        voxIdx(end+1) = sub2ind(scan_fov,tmpX(xi),tmpY(yi),tmpZ(zi)); %#ok<SAGROW>
                    end
                end
            end
        end
        
        % remove those indices which are not included in mask
        voxIdx(~ismember(voxIdx,M)) = []; %#ok<SAGROW>
               
        % convert whole brain indices to mask based indices
        sl_vox{vox} = find(ismember(M,voxIdx));
        
        % caluculate number of voxels in searchlight
        percentIncluded = numel(sl_vox{vox})/sum(reshape(sphere==1,1,[]));
        if percentIncluded < 0.6 % if less than 60% of max
            goodSL(vox) = false;
        end        
    end       
    
    % clean up
    clear percentIncluded voxIdx xi yi zi tmpX tmpY tmpZ vox

    % remove searchlights with less than threshold number of voxels
    sl_vox(goodSL==false) = [];
    
    % run LDt analysis
    RDM_ldt = rsa.stat.fisherDiscrTRDM_searchlight(X.a,Y.a,X.b,Y.b,1:8,sl_vox,model_rdm,false);

    % save raw output
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-rawRDM.mat'],'RDM_ldt')
    
    % clean up
    clear X Y sl_vox
    
    % add z-value to rdmBrain
    rdmBrain(M(goodSL)) = mean(cat(2,RDM_ldt.ats,RDM_ldt.bts),2);

    % load template nifti
    filename = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa/mask.nii'];        
    V = load_untouch_nii(filename);

    % change filename, datatype, and image
    V.fileprefix = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-searchlight'];
    V.hdr.dime.datatype = 64;
    V.img = rdmBrain;

    % save z-value brain
    save_untouch_nii(V,[V.fileprefix,'.nii']);  
    
    % update command line
    tElapse  = toc;
    tPerLoop = tElapse / subj;
    loopsRem = n_subj - subj;
    timeRem  = (tPerLoop * loopsRem)/3600;
    fprintf('\nSubject %1.0f of %1.0f complete...\nApproximate time remaining: %1.0f hours...\n',subj,n_subj,timeRem)
           
    % tidy
    clear V filename rdmBrain avgZ mN i goodSL M mask_idx RDM_ldt subjHandle
end

% clean up
clear subj tElapse tPerLoop loopsRem timeRem

%% Normalise and Smooth
% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % normalise functional
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb      = [-78 -112 -70; 78 76 85];    
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox     = [3 3 4];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp  = 4;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def         = {[dir_root,'bids_data/derivatives/',subj_handle,'/anat/y_',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample    = {[dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-searchlight.nii,1']};

    % smooth
    matlabbatch{2}.spm.spatial.smooth.fwhm                      = [8 8 8];
    matlabbatch{2}.spm.spatial.smooth.dtype                     = 0;
    matlabbatch{2}.spm.spatial.smooth.im                        = 1;
    matlabbatch{2}.spm.spatial.smooth.prefix                    = 's';
    matlabbatch{2}.spm.spatial.smooth.data                      = {[dir_root,'bids_data/derivatives/',subj_handle,'/rsa/w',subj_handle,'_task-rf_rsa-searchlight.nii,1']};
    
    % run batch
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
end

%% Run Second-Level Statistics
% predefine cell for searchlight image files
rMapFiles = cell(n_subj,2);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subj_handle = sprintf('sub-%02.0f',subj);
    
    % get searchlight images
    rMapFiles{subj,1}  = [dir_root,'bids_data/derivatives/',subj_handle,'/rsa/sw',subj_handle,'_task-rf_rsa-searchlight.nii,1'];
end

% create second-level glms
matlabbatch{1}.spm.stats.factorial_design.dir                       = {[dir_root,'bids_data/derivatives/group/rsa/visual-ers/']};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans              = rMapFiles(:,1);
matlabbatch{1}.spm.stats.factorial_design.cov                       = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov                 = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none        = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im                = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em                = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm           = 1;
spm_jobman('run',matlabbatch)
clear matlabbatch subjHandle subj

% estimate model
matlabbatch{1}.spm.stats.fmri_est.spmmat            = {[dir_root,'bids_data/derivatives/group/rsa/visual-ers/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals   = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1;

spm_jobman('run',matlabbatch)
clear matlabbatch

% define contrasts
matlabbatch{1}.spm.stats.con.spmmat(1)                  = {[dir_root,'bids_data/derivatives/group/rsa/visual-ers/SPM.mat']};   
matlabbatch{1}.spm.stats.con.delete                     = 0;    
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name       = 'within>between';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec     = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none';  

spm_jobman('run',matlabbatch)
clear matlabbatch

%% Extract Metrics for Visualisation
% combine visual cluster and save
combine_spm_cluster([dir_root,'bids_data/derivatives/group/rsa/visual-ers/'])

% load SPM details
load([dir_root,'bids_data/derivatives/group/rsa/visual-ers/SPM.mat'])

% extract subject values for each cluster
[betas,d] = extract_sample_points([dir_root,'bids_data/derivatives/group/rsa/visual-ers/'],SPM);

% save betas as table
tbl = array2table(betas','VariableNames',{'Occipital','Central'});
writetable(tbl,[dir_repos,'data/ers_betas.csv'],'Delimiter',',')

% save effect size as table
tbl = array2table(d','VariableNames',{'Occipital','Central'});
writetable(tbl,[dir_repos,'data/ers_cohensD.csv'],'Delimiter',',')
