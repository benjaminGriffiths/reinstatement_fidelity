function source = remove_source_electrodes(source,mask)

% get index of channels to include
i   = mask.inside(:);
j   = mask.insideWB(:);
idx = i(j==1);

% remove labels
source.label = source.label(idx);

% remove data from trials
for trl = 1 : numel(source.trial)
    source.trial{trl} = source.trial{trl}(idx,:);
end