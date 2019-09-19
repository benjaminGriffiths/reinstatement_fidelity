function data = recode_trlinfo(data)

% prepare matrix
trlinfo = nan(numel(data.trialinfo),6);

% cycle through each trial
for trl = 1 : numel(data.trialinfo)
    
    % extract details of interest
    trlinfo(trl,1) = data.trialinfo{trl}.recalled;
    trlinfo(trl,2) = data.trialinfo{trl}.confidence;
    trlinfo(trl,3) = data.trialinfo{trl}.trl_at_encoding;
    trlinfo(trl,4) = data.trialinfo{trl}.trl_at_retrieval;
    trlinfo(trl,5) = strcmpi(data.trialinfo{trl}.operation,'encoding');
    trlinfo(trl,6) = strcmpi(data.trialinfo{trl}.modality,'visual');
end

% repackage
data.trialinfo = trlinfo;
