function [freq] = get_rsa_timefrequency_correlation(data,operation,modality,si)

% define time frequency analysis parameters based on data probability
if numel(data.label) > 200
    fprintf('More than 200 channels detected, assuming data is in source space...\n')
    toi = -1:0.05:3;
    foi = 3:0.5:40;
else
    fprintf('Less than 200 channels detected, assuming data is channel level...\n')
    toi = -1:0.05:3;
    foi = 3:0.5:40;
end

% predefine conditional arrays to include all trials
operation_to_include = zeros(numel(data.trial),1);
modality_to_include  = zeros(numel(data.trial),1);

% cycle through each trial
for trl = 1 : numel(data.trial)

    % if encoding is requested#
    if ~isempty(operation)
        
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
    
end

% select data
cfg                 = [];
cfg.trials          = operation_to_include == 1 & modality_to_include == 1;
data                = ft_selectdata(cfg,data);

% predefine arrays for memory performance and trl number
mem_performance = zeros(numel(data.trial),1);
trl_no = zeros(numel(data.trial),1);

% cycle through each trial
for trl = 1 : numel(data.trial)

    % get memory
    mem_performance(trl) = data.trialinfo{trl}.recalled;
    
    % get trial number
    trl_no(trl) = data.trialinfo{trl}.(['trl_at_',operation]);    
end

% fix trial numbers
trl_no(trl_no>48&trl_no<=96)    = trl_no(trl_no>48&trl_no<=96) - 48;
trl_no(trl_no>96&trl_no<=144)   = trl_no(trl_no>96&trl_no<=144) - 48;
trl_no(trl_no>144)              = trl_no(trl_no>144) - 96;

% if retrieval, kick out forgotten trials
if strcmpi(operation,'retrieval')
    trl_no(mem_performance~=1)      = [];
end

% get similarity index for these trials
si = si(trl_no);

% define time-frequency parameters
cfg                 = [];
cfg.keeptrials      = 'yes';
cfg.method          = 'wavelet';
cfg.width           = 6;
cfg.output          = 'pow';
cfg.toi             = toi;
cfg.foi             = foi;
cfg.pad             = 'nextpow2';

% if retrieval, kick out forgotten trials
if strcmpi(operation,'retrieval')    
    cfg.trials = find(mem_performance==1);
end

% get time frequency
freq                = ft_freqanalysis(cfg, data);
clear data

% get time-averaged mean and standard deviation of power for each channel and frequency
raw_pow = mean(freq.powspctrm,4);
avg_pow = mean(raw_pow,1);
std_pow = std(raw_pow,[],1);

% replicate matrices to match freq.powspctrm
avg_pow = repmat(avg_pow,[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);
std_pow = repmat(std_pow,[size(freq.powspctrm,1) 1 1 size(freq.powspctrm,4)]);

% z-transform power
freq.powspctrm = (freq.powspctrm - avg_pow) ./ std_pow;
clear raw_pow avg_pow std_pow

% if encoding, get event-related change
if strcmpi(operation,'encoding')    
    
    % get pre- and post- power
    prepow  = mean(freq.powspctrm(:,:,:,freq.time<=0),4);
    postpow = mean(freq.powspctrm(:,:,:,freq.time>=0.5&freq.time<=1.5),4);
    
    % add difference to freq
    freq.powspctrm = postpow - prepow;
    
    % update time
    freq.time = 1;
end

% predefine matrix to hold correlation data
r = nan(size(freq.powspctrm,2),size(freq.powspctrm,3),size(freq.powspctrm,4));

% cycle through each channel
for chan = 1 : size(freq.powspctrm,2)

    % get parameters to correlate
    X = freq.powspctrm(:,chan,:,:);
    Y = repmat(si',[1 1 size(freq.powspctrm,3) size(freq.powspctrm,4)]);
    
    % correlate time-frequency data with similarity index
    numer = sum((X-mean(X,1)).*(Y-mean(Y,1)));
    denom = sqrt(sum((X-mean(X,1)).^2) .* sum((Y-mean(Y,1)).^2));
    r(chan,:,:) = atanh(squeeze(numer ./ denom));
end

% add correlation into freq structure
freq            = rmfield(freq,{'cumtapcnt','trialinfo'});
freq.powspctrm  = r;
freq.dimord     = 'chan_freq_time';
freq.cfg        = [];
