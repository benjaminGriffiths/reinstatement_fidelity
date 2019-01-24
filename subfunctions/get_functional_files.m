function ims = get_functional_files(subj_dir,prefix,copylocal)

if nargin ~= 3; copylocal = false; end

% get functional files
nii_files = dir([subj_dir,'func/',prefix,'*.nii']);

% if copy local is requested
if copylocal
    
    % make sure directory exists
    if ~exist('C:/tmp_mri/','file'); mkdir('C:/tmp_mri/'); end
    
    % cycle through each run
    for run = 1 : numel(nii_files)
        
        % copy files to temporary
        copyfile([subj_dir,'func/',nii_files(run).name],['C:/tmp_mri/',nii_files(run).name]);
        
        % get 4D scans in SPM-readable string format
        ims_per_run{run,1} = cellstr(spm_select('expand',['C:/tmp_mri/',nii_files(run).name]));
    end
    
    % get 4D scans in SPM-readable string format
    ims = {cat(1,ims_per_run{:})};
    
else
    % cycle through each run
    for run = 1 : numel(nii_files)

        % get 4D scans in SPM-readable string format
        ims_per_run{run,1} = cellstr(spm_select('expand',[subj_dir,'func/',nii_files(run).name]));
    end

    % get 4D scans in SPM-readable string format
    ims = {cat(1,ims_per_run{:})};    
end

