function freq = get_fooof(data,modality,operation)

% check variables
fprintf('organising data...\n')

% predefine conditional arrays to include all trials
modality_to_include   = zeros(numel(data.trial),1);
operation_to_include  = zeros(numel(data.trial),1);

% cycle through each trial
for trl = 1 : numel(data.trial)

    % mark trials that match specified operation
    modality_to_include(trl) = strcmpi(data.trialinfo{trl}.modality,modality);
    operation_to_include(trl) = strcmpi(data.trialinfo{trl}.operation,operation);
end

% get conditions which fall into both categories
comb_to_include = modality_to_include & operation_to_include;

% select trials of interest
data.trial      = data.trial(comb_to_include==1);
data.time       = data.time(comb_to_include==1);
data.trialinfo  = data.trialinfo(comb_to_include==1);
clear trl modality_to_include modality

% define key variables
ntrl    = numel(data.trial);
nchan   = numel(data.label);
ntime   = numel(data.time{1});

% extract signal
signal = permute(reshape(cell2mat(data.trial),[nchan,ntime,ntrl]),[2 1 3]);

% get time-limited power
postsig = signal(data.time{1}>=0.5&data.time{1}<=1.5,:,:);
presig  = signal(data.time{1}>=-1&data.time{1}<=0,:,:);
signal  = cat(3,postsig,presig);
clear postsig presig

% get PSD
fprintf('calculating PSD...\n');
[psd,freqs] = pwelch(signal(:,:),100,[],[],100);

% predefine output data
peak_amp   = nan(ntrl*2*nchan,1);
peak_freq  = nan(ntrl*2*nchan,1);
frac_slope = nan(ntrl*2*nchan,2);
auc        = nan(ntrl*2*nchan,3);

% cycle through each trl
fprintf('getting trial FOOOF...\n');
parfor i = 1 : size(psd,2)

    % run FOOOF
    modfit = fooof(freqs',psd(:,i)',[5 40],struct('peak_width_limits',[0.8 3],'max_n_peaks',8),true);

    % get peak amplitude/freq
    peaks   = modfit.peak_params(modfit.peak_params(:,1)<15,:);
    
    % if peak found
    if ~isempty(peaks)
        [~,idx]         = max(peaks(:,2));
        peak_amp(i,1)   = peaks(idx,1);
        peak_freq(i,1)  = peaks(idx,2);
    end
        
    % extract slope parameters
    frac_slope(i,:) = modfit.background_params;
    
    % get oscillatory spectrum
    osc_spc = modfit.fooofed_spectrum - modfit.bg_fit;
    
    % get area under slope
    auc(i,:) = [trapz(modfit.bg_fit),trapz(osc_spc(modfit.freqs>=8 & modfit.freqs<=14)),trapz(osc_spc(modfit.freqs>=18 & modfit.freqs<=30))];
end

% tidy workspace
clear freqs psd signal freq_idx

% reshape output
peak_amp    = reshape(peak_amp,[nchan,ntrl*2])';
peak_freq   = reshape(peak_freq,[nchan,ntrl*2])';
frac_slope  = permute(reshape(frac_slope,[nchan,ntrl*2,2]),[2 1 3]);
auc         = permute(reshape(auc,[nchan,ntrl*2,3]),[2 1 3]);

% reshape output
powspctrm   = cat(3,peak_amp,peak_freq,auc(:,:,2:3));
slope       = cat(3,frac_slope,auc(:,:,1),nan(size(peak_amp)));
clear peak_amp peak_freq frac_slope auc

% package data
freq = struct('time',1,...
              'freq',[1 2 10 20],...
              'label',{data.label},...
              'powspctrm',powspctrm,...
              'slope',slope,...
              'trialinfo',{repmat(data.trialinfo,[2 1])},...
              'dimord','rpt_chan_freq_time');
                   
% tidy workspace
clear peak_amp powspctrm slope nandx frac_slope nchan ntime ntrl data
fprintf('FOOOF done...\n')
