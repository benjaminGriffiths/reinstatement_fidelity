function freq = get_wavelet(data)

% get data dimensions
nchan = numel(data.label);
ntime = numel(data.time{1});
ntrl  = numel(data.trial);

% get trial types
trialtypes = false(ntrl,2);
for trl = 1 : ntrl
    trialtypes(trl,1) = strcmpi(data.trialinfo{trl}.modality,'visual');
    trialtypes(trl,2) = strcmpi(data.trialinfo{trl}.operation,'encoding');
end

% extract trial data
dat     = reshape(cell2mat(data.trial),[nchan ntime ntrl]);
time    = reshape(cell2mat(data.time),[1 ntime ntrl]);

% adjust freqs as neccesary
pad = 10.24;
foi = round((3:30) .* pad) ./ pad;
toi = -1.5:0.05:3;
  
% predefine wavelet outputs
spectrum = single(nan(ntrl,nchan,numel(foi),numel(toi)));
fprintf('running wavelet analysis...\n')

% cycle through channels and convolve
chkpoints = round(linspace(0,ntrl,11));
for trl = 1 : ntrl

    % get TFR
    spectrum(trl,:,:,:) = single(ft_specest_wavelet(dat(:,:,trl),time(:,:,trl),...
                                  'timeoi',toi,'pad',pad,'padtype','zero',...
                                  'width', 6, 'gwidth',3,'freqoi',foi,...
                                  'polyorder',0,'verbose',false));
                              
    % if trl is a checkpoint
    if any(ismember(chkpoints,trl))
        fprintf('%2.0f%% complete...\n',(find(ismember(chkpoints,trl))-1)*10)
    end
end
                                    
% convert complex result to power
spectrum = abs(spectrum).^2;

% baseline correct each trial type
for i = 1 : 2
    for j = 1 : 2
        
        % get idx
        idx = (trialtypes(:,1) == i-1) & (trialtypes(:,2) == j-1);
        
        % get time-averaged mean and standard deviation of power for each channel and frequency
        raw_pow = mean(spectrum(idx,:,:,:),4);
        avg_pow = mean(raw_pow,1);
        std_pow = std(raw_pow,[],1);

        % replicate matrices to match freq.powspctrm
        avg_pow = repmat(avg_pow,[sum(idx) 1 1 size(spectrum,4)]);
        std_pow = repmat(std_pow,[sum(idx) 1 1 size(spectrum,4)]);

        % z-transform power
        spectrum(idx,:,:,:) = (spectrum(idx,:,:,:) - avg_pow) ./ std_pow;
        clear raw_pow avg_pow std_pow
    end
end

% package data
fprintf('done...\n')
freq = struct('cfg',[],...
              'label',{data.label},...
              'freq',foi,...
              'time',toi,...
              'powspctrm',spectrum,...
              'trialinfo',{data.trialinfo},...
              'dimord','rpt_chan_freq_time');
