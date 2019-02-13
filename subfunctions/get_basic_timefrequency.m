function freq = get_basic_timefrequency(data,operation,modality,resolution)

% define time frequency analysis parameters based on data probability
if numel(data.label) > 200
    fprintf('More than 200 channels detected, assuming data is in source space...\n')
    toi = -1:0.05:3;
    foi = 8:30;
else
    fprintf('Less than 200 channels detected, assuming data is channel level...\n')
    toi = -1:0.05:3;
    foi = 3:0.5:40;
end

if nargin == 4
    toi = resolution.toi;
    foi = resolution.foi;
end

% predefine conditional arrays to include all trials
operation_to_include = zeros(numel(data.trial),1);
modality_to_include  = zeros(numel(data.trial),1);

% predefine arrays for memory performance
mem_performance = zeros(numel(data.trial),1);

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
    
    % get memory
    mem_performance(trl) = data.trialinfo{trl}.recalled;
end

% select data
cfg                 = [];
cfg.trials          = operation_to_include == 1 & modality_to_include == 1;
data                = ft_selectdata(cfg,data);

% select memory performance details
mem_performance = mem_performance(operation_to_include == 1 & modality_to_include == 1);

% get time-frequency parameters
cfg                 = [];
cfg.keeptrials      = 'yes';
cfg.method          = 'wavelet';
cfg.width           = 6;
cfg.output          = 'pow';
cfg.toi             = toi;
cfg.foi             = foi;
cfg.pad             = 'nextpow2';
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
        
    % resort trialinfo
    new_info = zeros(numel(freq.trialinfo),3);
    for i = 1 : numel(freq.trialinfo)
        new_info(i,:) = [freq.trialinfo{i}.trl_at_encoding mem_performance(i) freq.trialinfo{i}.confidence];
    end
    
elseif strcmpi(operation,'retrieval')       
    
    % get memory difference
    mem_diff = mean(mean(freq.powspctrm(mem_performance == 1,:,:)) - mean(freq.powspctrm(mem_performance == 0,:,:)),3);
    
    % kick out forgotten trials
    cfg         = [];
    cfg.latency = [0.5 1.5];
    freq        = ft_selectdata(cfg,freq);    
    
    % resort trialinfo
    new_info = zeros(numel(freq.trialinfo),3);
    for i = 1 : numel(freq.trialinfo)
        new_info(i,:) = [freq.trialinfo{i}.trl_at_retrieval mem_performance(i) freq.trialinfo{i}.confidence];
    end
        
    freq.X          = mem_diff;
end

% add correlation into freq structure
freq            = rmfield(freq,'cumtapcnt');
freq.trialinfo  = new_info;
freq.cfg        = [];
