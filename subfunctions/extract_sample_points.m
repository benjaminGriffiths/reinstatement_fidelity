function [betas,d] = extract_sample_points(data_dir,SPM)

% get number of clusters
n_files = numel(dir([data_dir,'\clus*.nii']));

% cycle through each cluster
for c = 1 : n_files

    % load in cluster
    nii = load_nii([data_dir,'clus',num2str(c),'.nii']);

    % get 3 x m representation of cluster
    roi = [];
    [roi(1,:),roi(2,:),roi(3,:)] = ind2sub(size(nii.img),find(nii.img>0));

    % get average of cluster
    betas(c,:) = nanmean(spm_get_data(SPM.xY.P,roi),2); 
end

% calculate cohens d for each cluster
for b = 1 :  n_files
    
    % calculate cohen's dz
    X = betas(b,:);
    numer = nanmean(X);
    denom = sqrt(sum((X - nanmean(X)).^2) ./ (numel(X)-1));
    d(b,1) = numer ./ denom;
end
