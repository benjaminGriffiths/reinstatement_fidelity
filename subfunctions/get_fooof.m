function [freq,spec,h] = get_fooof(data,modality)

% check variables
fprintf('organising data...\n')

% predefine conditional arrays to include all trials
modality_to_include  = zeros(numel(data.trial),1);

% cycle through each trial
for trl = 1 : numel(data.trial)

    % mark trials that match specified operation
    modality_to_include(trl) = strcmpi(data.trialinfo{trl}.modality,modality);
end

% select trials of interest
data.trial      = data.trial(modality_to_include==1);
data.time       = data.time(modality_to_include==1);
data.trialinfo  = data.trialinfo(modality_to_include==1);
clear trl modality_to_include modality

% define key variables
ntrl    = numel(data.trial);
nchan   = numel(data.label);
ntime   = numel(data.time{1});

% extract signal
signal = permute(reshape(cell2mat(data.trial),[nchan,ntime,ntrl]),[2 1 3]);

% get PSD
fprintf('calculating PSD...\n');
[psd,freqs] = pwelch(signal(:,:),100,[],[],100);
psd = reshape(psd,[size(psd,1),nchan,ntrl]);

% get trial averaged PSD
psd_avg = mean(psd,3);

% pre-define matrices
peak_freq   = nan(nchan,1);
flat_spec   = nan(nchan,floor(size(psd,1)/2));
freq_spec   = nan(nchan,floor(size(psd,1)/2));

% cycle through each channel
fprintf('getting average FOOOF...\n');
parfor i = 1 : nchan

    % run FOOOF
    modfit = fooof(freqs',psd_avg(:,i)',[5 30],struct('peak_width_limits',[1 4],'max_n_peaks',5),true);

    % find peaks
    peaks = modfit.peak_params;
    peaks(peaks(:,1)>=15 & peaks(:,1)<=17,:) = []; % cut TE artifact

    % find largest peak
    if ~isempty(peaks)
        [~,idx] = max(peaks(:,2));
        peak_freq(i,1) = peaks(idx,1);
    end

    % store flattened spectrum
    flat_spec(i,:) = modfit.fooofed_spectrum-modfit.bg_fit;   
    freq_spec(i,:) = modfit.freqs;            
end

% tidy workspace
clear psd psd_avg freqs

% get nan index and kick out trials
nandx = ~isnan(peak_freq);

% get peak frequency indices
freq_idx        = nan(sum(nandx),2);
freq_idx(:,1)   = knnsearch(freq_spec(1,:)',peak_freq(sum(nandx)))-3;
freq_idx(:,2)   = knnsearch(freq_spec(1,:)',peak_freq(sum(nandx)))+3;

% replicate peak frequency data
freq_idx = repmat(freq_idx,[ntrl*2,1]);

% normalise flattened spec.
spec.pow = (flat_spec - repmat(min(flat_spec,[],2),[1 size(flat_spec,2)])) ./ ...
            (repmat(max(flat_spec,[],2),[1 size(flat_spec,2)]) - repmat(min(flat_spec,[],2),[1 size(flat_spec,2)]));
spec.freq = freq_spec(1,:);
spec.peak = peak_freq;

% get time-limited power
postsig = signal(data.time{1}>=0.5&data.time{1}<=1.5,:,:);
presig  = signal(data.time{1}>=-1&data.time{1}<=0,:,:);
signal  = cat(3,postsig,presig);

% drop non-peak trials
signal = signal(:,nandx,:);

% tidy workspace
clear flat_spec freq_spec postsig presig peak_freq

% get PSD
fprintf('recalculating PSD...\n');
[psd,freqs] = pwelch(signal(:,:),100,[],[],100);

% predefine output data
peak_amp   = nan(ntrl*2*sum(nandx),1);
frac_slope = nan(ntrl*2*sum(nandx),1);

% cycle through each trl
fprintf('getting trial FOOOF...\n');
parfor i = 1 : size(psd,2)

    % run FOOOF
    modfit = fooof(freqs',psd(:,i)',[5 30],struct('peak_width_limits',[1 4],'max_n_peaks',5),true);

    % define peak index
    tmp_idx = freq_idx(i,:); % handle parallel broadcast issue
    if tmp_idx(1)<1; tmp_idx(1)=1; end
    if tmp_idx(2)>numel(modfit.freqs); tmp_idx(2)=numel(modfit.freqs); end
    
    % get peak amplitude
    flat_tmp = modfit.fooofed_spectrum-modfit.bg_fit;
    peak_amp(i,:) = mean(flat_tmp(tmp_idx(1):tmp_idx(2)));

    % get slope
    frac_slope(i) = modfit.background_params(2);
end

% tidy workspace
clear freqs psd signal freq_idx

% reshape output
peak_amp = reshape(peak_amp,[sum(nandx),ntrl*2]);
frac_slope = reshape(frac_slope,[sum(nandx),ntrl*2]);

% create new vars including nan trials
powspctrm   = nan(nchan,ntrl*2);
slope       = nan(nchan,ntrl*2);

% move data into nan-vars
powspctrm(nandx,:)  = peak_amp;
slope(nandx,:)      = frac_slope;

% package data
freq = struct('time',1,...
              'freq',8,...
              'label',{data.label},...
              'powspctrm',powspctrm',...
              'slope',slope',...
              'trialinfo',{repmat(data.trialinfo,[2 1])},...
              'dimord','rpt_chan_freq_time');
                   
          
% tidy workspace
clear peak_amp powspctrm slope nandx frac_slope nchan ntime ntrl
fprintf('FOOOF done...\n')

% plot results
if nargout == 3
    
    h=figure('position',[100 100 800 500]); hold on
    
    subplot(2,2,1); hold on
    shadedErrorBar(spec.freq,mean(spec.pow),sem(spec.pow));
    xlabel('Freq. (Hz)'); ylabel('Power (a.u.)'); title('Oscillatory Power')
    
    subplot(2,2,2); hold on
    histogram(spec.peak(:),100,'Normalization','probability','BinLimits',[5,30])
    xlabel('Freq. (Hz)'); ylabel('Probability'); title('Bin Peak')
    
    subplot(2,2,3); hold on
    histogram(freq.powspctrm(:),100,'Normalization','probability')
    xlabel('Power (a.u.)'); ylabel('Probability'); title('Peak Power');
    
    subplot(2,2,4); hold on
    histogram(freq.slope(:),100,'Normalization','probability')
    xlabel('Slope Coefficent'); ylabel('Probability'); title('Peak Power');
end


