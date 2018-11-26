function ims_out = update_functional_files(ims,prefix)

% cycle through each file
for i = 1 : numel(ims{1})
    
    % get length of filename
    lenFile = numel(ims{1}{i}(1:end-2));
    
    % break down file
    [path,name,ext] = fileparts(ims{1}{i}(1:lenFile));
    
    ims_out{i,1} = [path,'\',prefix,name,ext,ims{1}{i}(lenFile+1:end)];
end
