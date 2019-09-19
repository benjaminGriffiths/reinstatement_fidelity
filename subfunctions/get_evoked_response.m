function freq = get_evoked_response(data,operation,modality)

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
cfg.latency         = [-0.25 0.25];
data                = ft_selectdata(cfg,data);

% get time-frequency parameters
cfg                 = [];
cfg.keeptrials      = 'yes';
tml                 = ft_timelockanalysis(cfg,data);

% find peak across trials
avg = squeeze(mean(tml.trial,1));
base = mean(avg(:,tml.time<0),2);
avg = avg - repmat(base,[1 size(avg,2)]);
for chan = 1 : size(avg,1)
    [val,tmp] = findpeaks(avg(chan,tml.time>0));
    if numel(val) > 1
        [~,midx] = max(val);
        idx(chan,1) = tmp(midx);
    else
        idx(chan,1) = tmp;
    end
end

% extract peak amp. for each trial
tmp = tml.trial - repmat(base',[size(tml.trial,1) 1 size(tml.trial,3)]);
tmp = tmp(:,:,tml.time>0);
peak_amp = zeros(size(tml.trial,1),size(tml.trial,2));
for chan = 1 : size(tmp,2)
    peak_amp(:,chan) = tmp(:,chan,idx(chan));
end

% add correlation into freq structure
freq            = [];
freq.freq       = 1;
freq.time       = 1;
freq.label      = data.label;
freq.dimord     = 'rpt_chan_freq_time';
freq.powspctrm  = peak_amp;
freq.cfg        = [];
