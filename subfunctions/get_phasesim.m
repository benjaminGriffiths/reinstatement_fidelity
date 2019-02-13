function [freq,ntrls] = get_phasesim(data,operation,modality)

% define time/frequency resolution
toi = 0:0.02:1;
foi = 8;

% predefine conditional arrays to include all trials
operation_to_include = zeros(numel(data.trial),1);
modality_to_include  = zeros(numel(data.trial),1);

% predefine arrays for memory performance
mem_performance = zeros(numel(data.trial),1);

% cycle through each trial
for trl = 1 : numel(data.trial)

    % if encoding is requested#
    if strcmpi(operation,'encoding')
        
        % mark trials that do match specified operation
        operation_to_include(trl) = strcmpi(data.trialinfo{trl}.operation,operation);
     
    % otherwise include
    else; operation_to_include(trl) = 1;
    end
    
    % if encoding is requested
    if ~isempty(modality)
        
        % mark trials that do match specified operation
        modality_to_include(trl) = strcmpi(data.trialinfo{trl}.modality,modality);
     
    % otherwise include
    else; modality_to_include(trl) = 1;
    end
    
    % get memory
    mem_performance(trl) = data.trialinfo{trl}.recalled;
end

% select data
cfg                 = [];
cfg.trials          = operation_to_include == 1 & modality_to_include == 1;
data                = ft_selectdata(cfg,data);

% get time-frequency
cfg                 = [];
cfg.keeptrials      = 'yes';
cfg.method          = 'wavelet';
cfg.width           = 6;
cfg.output          = 'fourier';
cfg.toi             = toi;
cfg.foi             = foi;
cfg.pad             = 'nextpow2';
freq                = ft_freqanalysis(cfg, data);
clear data

% extract phase
phase = angle(freq.fourierspctrm);

% extract stim values
stimval = nan(numel(freq.trialinfo),1);
for trl = 1 : numel(freq.trialinfo)
    stimval(trl,1) = freq.trialinfo{trl}.stimulus_value;
end

% predefine observed and predicted model
obs_model   = [];
hyp_model   = [];
count       = 1;

% cycle through every trial
for i = 1 : size(phase,1)
    for j = i+1 : size(phase,1)
        
        % extract signal
        X = squeeze(phase(i,:,:,:));
        Y = squeeze(phase(j,:,:,:));
        
        % get single-trial phase-locking 
        obs_model(count,:) = abs(mean(exp(1i.*(X-Y)),2)); 
        
        % get prediction
        hyp_model(count,1) = stimval(i)==stimval(j);
        count = count + 1;        
    end
end

% predefine correlation vector
r = nan(size(obs_model,2),1);

% cycle through each channel
for chan = 1 : size(obs_model,2)
    
    % correlate observed and predicted model
    r(chan,1) = atanh(corr(obs_model(:,chan),hyp_model));
end

% get difference between means
freq            = rmfield(freq,{'fourierspctrm','cumtapcnt','cfg','trialinfo'});
freq.powspctrm  = r;
freq.dimord     = 'chan_freq_time';
freq.time       = 1;
freq.freq       = 1;
