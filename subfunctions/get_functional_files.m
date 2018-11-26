function ims = get_functional_files(subj_dir,prefix)

% get functional files
nii_files = dir([subj_dir,'func\',prefix,'*.nii']);

% cycle through each run
for run = 1 : numel(nii_files)
    
    % get 4D scans in SPM-readable string format
    ims_per_run{run,1} = cellstr(spm_select('expand',[subj_dir,'func\',nii_files(run).name]));
end

% get 4D scans in SPM-readable string format
ims = {cat(1,ims_per_run{:})};