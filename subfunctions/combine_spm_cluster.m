function combine_spm_cluster(cluster_dir)

% get cluster files
files = dir([cluster_dir,'clus*.nii']);

% load first file
nii = load_nii([files(1).folder,'\',files(1).name]);

% extract image
img = nii.img;

% cycle through remaining images
for i = 2 : numel(files)
    
    % load file
    nii = load_nii([files(i).folder,'\',files(i).name]);

    % add to image
    img = img + nii.img;
end

% create new nifti variable
nii_grand               = nii;
nii_grand.img           = uint8(img > 0);
nii_grand.fileprefix    = [cluster_dir,'\grand_cluster'];

% save grand cluster
save_nii(nii_grand,[nii_grand.fileprefix,'.nii'])
fprintf('Mask saved...\n')

% dilate cluster
dilate_mask([nii_grand.fileprefix,'.nii'],10,[3 3 4])

% load t image
nii = load_nii([cluster_dir,'\spmT_0001.nii']);

% mask t image
nii.img = nii.img .* single(nii_grand.img);

% rename t image
nii.fileprefix = [nii.fileprefix,'_masked'];

% save t image
save_nii(nii,[nii.fileprefix,'.nii'])
fprintf('Masked T saved...\n')

% reslice (to smooth)
reslice_nii([nii.fileprefix,'.nii'],[nii.fileprefix,'_smooth.nii'],[1 1 1]);
fprintf('Smoothed masked T saved...\n')