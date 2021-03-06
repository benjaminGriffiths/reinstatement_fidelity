function [freq,ntrls] = get_memdiff_timefrequency(data,operation,modality,resolution)

% define time frequency analysis parameters based on data probability
if numel(data.label) > 200
    fprintf('More than 200 channels detected, assuming data is in source space...\n')
    toi = 0.5:0.1:1.5;
    foi = 8:2:30;
else
    fprintf('Less than 200 channels detected, assuming data is channel level...\n')
    toi = -1:0.05:3;
    foi = 3:0.5:40;
end

% define custom resolution if specified
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

% get time-frequency
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

% split data by memory performance
cfg         = [];
cfg.trials  = find(mem_performance(operation_to_include == 1 & modality_to_include == 1) == 1);
mem1        = ft_selectdata(cfg,freq);

cfg         = [];
cfg.trials  = find(mem_performance(operation_to_include == 1 & modality_to_include == 1) == 0);
mem2        = ft_selectdata(cfg,freq);
clear freq

% get condition means
mem1        = ft_freqdescriptives([],mem1);
mem2        = ft_freqdescriptives([],mem2);

% get difference between means
cfg             = [];
cfg.parameter   = 'powspctrm';
cfg.operation   = 'subtract';
freq            = ft_math(cfg,mem1,mem2);
clear mem1 mem2

% count number of trials in each condition
ntrls = [sum(mem_performance(operation_to_include == 1 & modality_to_include == 1) == 1),...
         sum(mem_performance(operation_to_include == 1 & modality_to_include == 1) == 0)];
