%% RSA
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
scan_fov    = [64 64 32];                                   % scan field of view
scan_vox    = [3 3 4];                                      % scan voxel size
scan_search = 12;                                           % searchlight radius
scan_func   = {'_3_1','_4_1','_5_1','_6_1',...
               '_8_1','_9_1','_10_1','_11_1'};              % functional scan suffix

% add subfunctions
addpath([dir_repos,'subfunctions'])

%% Prepare Masks
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/',subj_handle,'/'];
    
    % make directory for masks
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/masks/'])
    
    % prepare deformation batch
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def        = {[dir_root,'bids_data/derivatives/',subj_handle,'/anat/iy_',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.util.defs.comp{1}.inv.space              = {[dir_root,'bids_data/',subj_handle,'/anat/',subj_handle,'_T1w.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.fnames             = {[dir_root,'bids_data/sourcedata/masks/whole_brain.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.weight             = {''};
    matlabbatch{1}.spm.util.defs.out{1}.push.savedir.saveusr    = {[dir_root,'bids_data/derivatives/',subj_handle,'/masks/']};
    matlabbatch{1}.spm.util.defs.out{1}.push.fov.file           = {[dir_root,'bids_data/derivatives/',subj_handle,'/func/meanua',subj_handle,'_task-rf_run-1_bold.nii']};
    matlabbatch{1}.spm.util.defs.out{1}.push.preserve           = 0;
    matlabbatch{1}.spm.util.defs.out{1}.push.fwhm               = [0 0 0];
    matlabbatch{1}.spm.util.defs.out{1}.push.prefix             = '';
    
    % run
    spm_jobman('run',matlabbatch)
    clear matlabbatch
    
    % move file to subject directory
    movefile([dir_root,'bids_data/derivatives/',subj_handle,'/masks/wwhole_brain.nii'],...
        [dir_root,'bids_data/derivatives/',subj_handle,'/masks/whole_brain.nii'])
    
    clear matlabbatch
end

%% Read Data
% cycle through each subject
for subj = 1 : n_subj
    
    % update command line
    fprintf('/n--- Working on Subject %d ---------/n',subj)
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/',subj_handle,'/'];
    
    % predefine matrix for mask data
    maskImg = zeros(1, prod(scan_fov));
    
    % load mask and add to matrix
    nii = load_untouch_nii([dir_root,'bids_data/derivatives/',subj_handle,'/masks/whole_brain.nii']);
    maskImg(1,:) = reshape(nii.img,1,[]);

    % predefine matrix for functional data
    scanVec = zeros(n_volumes.*numel(scan_func),numel(nii.img));
    
    % start scan counter
    scanCount = 1;
    
    % cycle through each run
    for i = 1 : numel(scan_func)
        
        % define functional filenames
        filename = [dir_root,'bids_data/derivatives/',subj_handle,'/func/',...
            'ua',subj_handle,'_task-rf_run-',num2str(i),'_bold.nii'];

        % read in nifti file
        nii = load_untouch_nii(filename);

        % cycle through each scan
        for j = 1 : n_volumes
            
            % extract image
            scanVec(scanCount,:) = reshape(nii.img(:,:,:,j),1,[]);
            
            % count scan as read in
            scanCount = scanCount + 1;            
        end
        
        % update command line
        fprintf('Run %d of %d read in.../n',i,numel(scan_func))
    end
    
    % clear up
    clear i j filename scanCount nii
    
    % save
    fprintf('/nSaving full volume.../n')
    mkdir([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/'])
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-fullVolume.mat'],'scanVec')
    
    % apply mask to scans
    fprintf('Masking data.../n')
    
    % define patterns for mask
    patterns = scanVec;

    % mask functional data (set all zero-element voxels in mask to zero)
    patterns(:,maskImg(1,:)==0) = 0;

    % locate "dead voxels" (functional data within mask that has zero values; e.g. where mask captures skull)
    deadIdx = any(patterns == 0);

    % set these voxels to zero across all scans
    patterns(:,deadIdx) = 0;
    
    % get a boolean vector of non-zero voxels in patterns matrix
    mask_idx = all(patterns~=0);
    
    % remove zero elements from patterns
    patterns(:,any(patterns == 0)) = [];
    
    % save
    fprintf('Saving masked volumes.../n')
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-maskedVolume.mat'],'patterns')
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-mask.mat'],'mask_idx')
    
    % clear excess variables
    clear scanVec nii deadIdx patterns mask_idx subjHandle maskImg
end

%% Mean Pattern Subtraction
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/',subj_handle,'/'];
    
    % load pattern data
    load([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-maskedVolume.mat'])
    
    % get mean pattern across all trials and replicate matrix
    meanPattern = repmat(mean(patterns,1),[size(patterns,1), 1]);

    % subtract mean pattern from data
    patterns = patterns - meanPattern;
    
    % save patterns
    save([dir_root,'bids_data/derivatives/',subj_handle,'/rsa/',subj_handle,'_task-rf_rsa-maskedDemeanedVolume.mat'],'patterns')
    
    % clean up
    clear subjHandle meanPattern patterns
end

%% Prepare GLMs and Data for Searchlight Analysis
% cycle through each subject
for subj = 1 : n_subj
    
    % define key subject strings
    subj_handle = sprintf('sub-%02.0f',subj);
    dir_subj = [dir_root,'bids_data/',subj_handle,'/'];
    
    % load pattern data
    load([dir_root,'data/fmri/rsa/data/',subjHandle,'/sl_volume_demeaned.mat'])
    
    % load SPM.mat
    load([dir_root,'data/fmri/rsa/data/',subjHandle,'/SPM.mat'])
    
    % load stimulus details
    load([dir_root,'data/fmri/rsa/data/',subjHandle,'/stim.mat'])
    
    % get design matrix (X) and split into two groups (Xa and Xb)
    X.raw = SPM.xX.X;
    
    % define conditions to calculate similarity between
    conditions = {'eBIKE','eFARM','eUNDERWATER','eWATERMILL','eACCORDIAN','eGUITAR','ePIANO','eTRUMPET',...
        'rBIKE','rFARM','rUNDERWATER','rWATERMILL','rACCORDIAN','rGUITAR','rPIANO','rTRUMPET'};
    
    % get matrix of trials by stimulus content
    trl_by_con = zeros(size(stim,1) ./ numel(conditions),numel(conditions));
    for i = 1 : numel(conditions)
        trl_by_con(:,i) = find(ismember(stim(:,2),conditions(i)));
    end
    
    % get GLM for nuisance regressors
    X.n = X.raw(:,max(trl_by_con(:))+1:end);
    
    % get condition average GLM
    for i = 1 : numel(conditions)
        X.c(:,i) = sum(X.raw(:,trl_by_con(:,i)),2);
    end
    
    % clean up
    clear conditions trl_by_con SPM stim
    
    % split GLM into two groups (train and test data) [excluding nuisance regressors]
    X.at = X.raw(1:size(X.raw)/2,1:n_trials);
    X.bt = X.raw((size(X.raw)/2)+1:end,n_trials+1:n_trials*2);
    
    % repeat for condition-average GLM
    X.aa = X.c(1:size(X.raw)/2,:);
    X.ba = X.c((size(X.raw)/2)+1:end,:);
    
    % add nuisance regressors (skipping constants)
    X.at(:,end+1:end+size(X.n,2)-1) = X.n(1:size(X.raw)/2,1:end-1);
    X.bt(:,end+1:end+size(X.n,2)-1) = X.n((size(X.raw)/2)+1:end,1:end-1);
    
    % repeat for condition average
    X.aa(:,end+1:end+size(X.n,2)-1) = X.n(1:size(X.raw)/2,1:end-1);
    X.ba(:,end+1:end+size(X.n,2)-1) = X.n((size(X.raw)/2)+1:end,1:end-1);
    
    % find zero-value rows
    X.a_badRow = all(X.at(:,1:n_trials)==0,2);
    X.b_badRow = all(X.bt(:,1:n_trials)==0,2);
    
    % remove zero-value rows
    X.at(X.a_badRow,:) = [];
    X.aa(X.a_badRow,:) = [];
    X.bt(X.b_badRow,:) = [];
    X.ba(X.b_badRow,:) = [];
    
    % cycle through A/B avg/trl combinations
    fn = {'at','aa','bt','ba'};
    for i = 1 : 4
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
    X.at(:,X.at_badCol) = [];
    X.aa(:,X.aa_badCol) = [];
    X.bt(:,X.bt_badCol) = [];
    X.ba(:,X.ba_badCol) = [];
    
    % get activation matrix (Y) and split into two groups (Ya and Yb)
    Y.raw = patterns;
    Y.a = Y.raw(1:size(Y.raw,1)/2,:);
    Y.b = Y.raw((size(Y.raw,1)/2)+1:end,:);

    % remove zero-value rows
    Y.a(X.a_badRow,:) = [];
    Y.b(X.b_badRow,:) = [];

    % save data
    save([dir_root,'data/fmri/rsa/data/',subjHandle,'/sl_volume_formatted.mat'],'X','Y')
    
    clear X Y i j fn patterns subjHandle
    
    fprintf('Subject %02.0f of %02.0f prepared.../n',subj,n_subj)
end

%% Define Searchlight Characteristics
% get the radius of searchlight in voxels
voxRadInSearchLight = scan_search ./ scan_vox;

% define distance from searchlight centre to perimeter in voxels
dist2Perimeter = floor(voxRadInSearchLight);

% create boolean searchlight sphere
[x,y,z] = meshgrid(-dist2Perimeter(1):dist2Perimeter(1),-dist2Perimeter(2):dist2Perimeter(2),-dist2Perimeter(3):dist2Perimeter(3));
sphere  = ((x*scan_vox(1)).^2+(y*scan_vox(2)).^2+(z*scan_vox(3)).^2)<=(scan_search^2);

% predefine model rdm
model_rdm{1} = nan(16,16);
model_rdm{2} = nan(16,16);

% set items that belong to the same category to -1
model_rdm{1}(9:12,1:4) = 1;
model_rdm{2}(13:16,5:8) = 1;

% set items that belong to a different category to 1
model_rdm{1}(9,1)  = -1;
model_rdm{1}(10,2) = -1;
model_rdm{1}(11,3) = -1;
model_rdm{1}(12,4) = -1;
model_rdm{2}(13,5) = -1;
model_rdm{2}(14,6) = -1;
model_rdm{2}(15,7) = -1;
model_rdm{2}(16,8) = -1;

% set diagonal to NaN
model_rdm{1}(logical(eye(size(model_rdm{1},1)))) = 0; 
model_rdm{2}(logical(eye(size(model_rdm{2},1)))) = 0; 

% clean up
clear x y z voxRadInSearchLight

%% Run Searchlight LDt Analysis
tic

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subjHandle = sprintf('subj%02.0f',subj);
    
    % load data
    load([dir_root,'data/fmri/rsa/data/',subjHandle,'/sl_volume_formatted.mat'])    
    load([dir_root,'data/fmri/rsa/data/',subjHandle,'/sl_mask_indexed.mat'])    
       
    % define empty brain map
    rdmBrain = zeros(scan_fov);
    
    % predefine voxels to be analysed
    sl_vox = cell(sum(mask_idx),1);
    
    % predefine vector of searchlights to be used
    goodSL = true(size(sl_vox));
    
    % get index of all voxels within mask
    M = find(mask_idx);
    
    % update command line
    fprintf('Defining searchlights.../n')
    
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
        if percentIncluded < 0.9 % if less than 90% of max
            goodSL(vox) = false;
        end        
    end       
    
    % clean up
    clear percentIncluded voxIdx xi yi zi tmpX tmpY tmpZ vox

    % remove searchlights with less than threshold number of voxels
    sl_vox(goodSL==false) = [];
    
    % define conditions to calculate similarity between
    conditions = {'eBIKE','eFARM','eUNDERWATER','eWATERMILL','eACCORDIAN','eGUITAR','ePIANO','eTRUMPET',...
        'rBIKE','rFARM','rUNDERWATER','rWATERMILL','rACCORDIAN','rGUITAR','rPIANO','rTRUMPET'};
    
    % run LDt analysis
    RDM_ldt = rsa.stat.fisherDiscrTRDM_searchlight(X.aa,Y.a,X.ba,Y.b,1:numel(conditions),sl_vox,model_rdm);

    % clean up
    clear X Y conditions model_rdm sl_vox
    
    % model names
    mN = {'visual','auditory'};
    
    for i = 1 : size(RDM_ldt.ats,2)
        % get mean corrcoef
        avgZ = mean(cat(2,RDM_ldt.ats(:,i),RDM_ldt.bts(:,i)),2);

        % add z-value to rdmBrain
        rdmBrain(M(goodSL)) = avgZ;

        % load template nifti
        filename = [dir_root,'data/fmri/preprocessing/subj',sprintf('%02.0f',subj),...
                '/uasubj',num2str(subj),scan_func{i},'_',sprintf('%05.0f',1),'.nii'];        
        V = load_untouch_nii(filename);

        % change filename, datatype, and image
        V.fileprefix = [dir_root,'data/fmri/rsa/data/',subjHandle,'/sl_ldt_',mN{i}];
        V.hdr.dime.datatype = 64;
        V.img = rdmBrain;

        % save z-value brain
        save_untouch_nii(V,[V.fileprefix,'.nii']);       
    end
    
    % update command line
    tElapse  = toc;
    tPerLoop = tElapse / subj;
    loopsRem = numel(n_subj) - subj;
    timeRem  = (tPerLoop * loopsRem)*3600;
    fprintf('/nSubject %1.0f of %1.0f complete.../nApproximate time remaining: %1.0f hours.../n',subj,n_subj,timeRem)
           
    % tidy
    clear V filename rdmBrain avgZ mN i goodSL M mask_idx RDM_ldt subjHandle
end

% clean up
clear subj tElapse tPerLoop loopsRem timeRem

%% Normalise and Smooth
% define images to analyse
img = {'sl_ldt_auditory.nii','sl_ldt_visual.nii'};

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subjHandle = sprintf('subj%02.0f',subj);
    
    % normalise functional
    matlabbatch{1}.spm.spatial.normalise.write.woptions.bb      = [-78 -112 -70; 78 76 85];    
    matlabbatch{1}.spm.spatial.normalise.write.woptions.vox     = [3 3 4];
    matlabbatch{1}.spm.spatial.normalise.write.woptions.interp  = 4;
    matlabbatch{1}.spm.spatial.normalise.write.subj.def         = {[dir_root,'data/fmri/preprocessing/',subjHandle,'/y_subj',num2str(subj),'_7_1.nii']};
    matlabbatch{1}.spm.spatial.normalise.write.subj.resample    = {[dir_root,'data/fmri/rsa/data/',subjHandle,'/',img{1},',1'];
                                                                   [dir_root,'data/fmri/rsa/data/',subjHandle,'/',img{2},',1']};

    % run batch
    spm_jobman('run',matlabbatch)
    
    % smooth
    matlabbatch{2}.spm.spatial.smooth.fwhm                      = [8 8 8];
    matlabbatch{2}.spm.spatial.smooth.dtype                     = 0;
    matlabbatch{2}.spm.spatial.smooth.im                        = 0;
    matlabbatch{2}.spm.spatial.smooth.prefix                    = 's';
    matlabbatch{2}.spm.spatial.smooth.data                      = {[dir_root,'data/fmri/rsa/data/',subjHandle,'/w',img{1},',1'];
                                                                   [dir_root,'data/fmri/rsa/data/',subjHandle,'/w',img{2},',1']};
    
    % run batch
    spm_jobman('run',matlabbatch)
    clear matlabbatch img subjHandle
end

%% Run Second-Level Statistics
% predefine cell for searchlight image files
rMapFiles = cell(n_subj,2);

% cycle through each subject
for subj = 1 : n_subj
    
    % define subject name
    subjHandle = sprintf('subj%02.0f',subj);
    
    % get searchlight images
    rMapFiles{subj,1}  = [dir_root,'data/fmri/rsa/data/',subjHandle,'/swsl_ldt_visual.nii'];
    rMapFiles{subj,2}    = [dir_root,'data/fmri/rsa/data/',subjHandle,'/swsl_ldt_auditory.nii'];
end

% create second-level glm
matlabbatch{1}.spm.stats.factorial_design.dir                       = {[dir_root,'data/fmri/rsa/stats/sl_visual']};
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
matlabbatch{1}.spm.stats.fmri_est.spmmat            = {[dir_root,'data/fmri/rsa/stats/sl_visual/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals   = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1;

spm_jobman('run',matlabbatch)
clear matlabbatch

% define contrasts
matlabbatch{1}.spm.stats.con.spmmat(1)                  = {[dir_root,'data/fmri/rsa/stats/sl_visual/SPM.mat']};   
matlabbatch{1}.spm.stats.con.delete                     = 0;    
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name       = 'VideoSimiliarity';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec     = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none';  

spm_jobman('run',matlabbatch)
clear matlabbatch

% create second-level glm
matlabbatch{1}.spm.stats.factorial_design.dir                       = {[dir_root,'data/fmri/rsa/stats/sl_auditory']};
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans              = rMapFiles(:,2);
matlabbatch{1}.spm.stats.factorial_design.cov                       = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov                 = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none        = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im                = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em                = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit            = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no    = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm           = 1;
spm_jobman('run',matlabbatch)
clear matlabbatch

% estimate model
matlabbatch{1}.spm.stats.fmri_est.spmmat            = {[dir_root,'data/fmri/rsa/stats/sl_auditory/SPM.mat']};
matlabbatch{1}.spm.stats.fmri_est.write_residuals   = 0;
matlabbatch{1}.spm.stats.fmri_est.method.Classical  = 1;

spm_jobman('run',matlabbatch)
clear matlabbatch

% define contrasts
matlabbatch{1}.spm.stats.con.spmmat(1)                  = {[dir_root,'data/fmri/rsa/stats/sl_auditory/SPM.mat']};   
matlabbatch{1}.spm.stats.con.delete                     = 0;    
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name       = 'AudioSimiliarity';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec     = 1;
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep    = 'none';  

spm_jobman('run',matlabbatch)
clear matlabbatch rMapFiles

%% Extract Data from Cluster
% load SPM details
load([dir_root,'data/fmri/rsa/stats/sl_visual/SPM.mat'])

% define cluster name
cluster_name = {'left_fusiform_181120','right_fusiform_181120','right_pof_181120'};

% cycle through each cluster
for c = 1 : 3

    % load in cluster
    nii = load_nii([dir_root,'data/fmri/rsa/stats/sl_visual/',cluster_name{c},'.nii']);

    % get 3 x m representation of cluster
    roi = [];
    [roi(1,:),roi(2,:),roi(3,:)] = ind2sub(size(nii.img),find(nii.img>0));

    % get average of cluster
    betas(c,:) = mean(spm_get_data(SPM.xY.P,roi),2); 
end

% rename each beta and add to table
leftFusiform    = betas(1,:)';
rightFusiform   = betas(2,:)';
rightPOF        = betas(3,:)';

% add to table
data = table(leftFusiform,rightFusiform,rightPOF);

% calculate cohens d
for b = 1 : 3
    
    % calculate cohen's dz
    X = betas(b,:);
    numer = mean(X);
    denom = sqrt(sum((X - mean(X)).^2) ./ (numel(X)-1));
    d(b,1) = numer ./ denom;
end

% export
writetable(data,[dir_repos,'data/searchlight_betas.csv'],'Delimiter',',')

