function [audio_bool,hit_bool] = get_splits_eeg(freq)

% get splits in data
hit_bool   = nan(size(freq.trialinfo));
audio_bool = nan(size(freq.trialinfo));
for trl = 1 : numel(freq.trialinfo)

    % identify whether the trial involved an auditory stimulus
    if freq.trialinfo{trl}.dynStimType == 1                  
        audio_bool(trl) = 1;
    else
        audio_bool(trl) = 0;
    end

    % identify whether the stimulus was correctly recalled with confidence
    if freq.trialinfo{trl}.hitBoolean == 1 && freq.trialinfo{trl}.confidenceScore > 1
        hit_bool(trl) = 1;            

    elseif freq.trialinfo{trl}.hitBoolean == 0 || freq.trialinfo{trl}.confidenceScore == 1
        hit_bool(trl) = 0;

    end
end