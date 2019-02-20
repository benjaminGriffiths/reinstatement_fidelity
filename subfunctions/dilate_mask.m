function dilate_mask(filename,searchlight,voxelsize)

% load mask
nii = load_nii(filename);

% create copy of image
img = zeros(size(nii.img));

% get the radius of searchlight in voxels
voxRadInSearchLight = searchlight ./ voxelsize;

% define distance from searchlight centre to perimeter in voxels
dist2Perimeter = ceil(voxRadInSearchLight);

% create boolean searchlight sphere
[x,y,z] = meshgrid(-dist2Perimeter(1):dist2Perimeter(1),-dist2Perimeter(2):dist2Perimeter(2),-dist2Perimeter(3):dist2Perimeter(3));
sphere  = ((x*voxelsize(1)).^2+(y*voxelsize(2)).^2+(z*voxelsize(3)).^2)<=(searchlight^2);

% clean up
clear x y z voxRadInSearchLight

% define field of view
fov = size(nii.img);

% cycle through every voxel
for vox = 1 : numel(img)

    % move on if voxel is not in mask
    if nii.img(vox)<1; continue; end
    
    % get subscript co-ordinates of searchlight centre
    [x,y,z] = ind2sub(size(nii.img),vox);

    % define box which houses spherical searchlight
    tmpX = x-dist2Perimeter(1):x+dist2Perimeter(1);
    tmpY = y-dist2Perimeter(2):y+dist2Perimeter(2);
    tmpZ = z-dist2Perimeter(3):z+dist2Perimeter(3);

    % remove co-ordinates less than or equal to zero, and greater than FOV
    tmpX(tmpX<=0) = []; tmpX(tmpX>fov(1)) = [];
    tmpY(tmpY<=0) = []; tmpY(tmpY>fov(2)) = [];
    tmpZ(tmpZ<=0) = []; tmpZ(tmpZ>fov(3)) = [];

    % cycle through 3D space
    for xi = 1:numel(tmpX)
        for yi = 1:numel(tmpY)
            for zi = 1:numel(tmpZ)       
                if sphere(xi,yi,zi)       
                    % mark voxel in new image
                    img(tmpX(xi),tmpY(yi),tmpZ(zi)) = 1;   
                end
            end
        end
    end
end 

% save new mask
nii.img         = uint8(img);
nii.fileprefix  = [nii.fileprefix,'_dilated'];
save_nii(nii,[nii.fileprefix,'.nii'])
