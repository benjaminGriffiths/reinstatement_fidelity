function freq = get_phasesim(data,operation,modality)

% define time/frequency resolution
if strcmpi(operation,'encoding')
    toi = data.time{1}(data.time{1}>=0&data.time{1}<=1);
    foi = 8;
else    
    toi = data.time{1}(data.time{1}>=0&data.time{1}<=2);
    foi = 8;
end

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

% extract stim and operation values
stimval = nan(numel(freq.trialinfo),1);
opval = nan(numel(freq.trialinfo),1);
for trl = 1 : numel(freq.trialinfo)
    stimval(trl,1) = freq.trialinfo{trl}.stimulus_value;
    opval(trl,1) = strcmpi(freq.trialinfo{trl}.operation,'encoding');
end

% predefine observed and predicted model
obs_model   = [];
hyp_model   = [];
count       = 1;

% change method based on whether encoding-encoding similarity or
% encoding-retreival similarity
if strcmpi(operation,'encoding')

    % cycle through every trial
    for i = 1 : size(phase,1)
        for j = i+1 : size(phase,1)

            % extract signal
            X = squeeze(phase(i,:,:,:));
            Y = squeeze(phase(j,:,:,:));

            % get single-trial phase-locking 
            obs_model(count,:) = mean(cos(X).*cos(Y) + sin(X).*sin(Y),2);   

            % get prediction
            hyp_model(count,1) = stimval(i)==stimval(j);
            count = count + 1;        
        end
    end
    
else
    
    % split encoding and retrieval data
    enc_phase = phase(opval==1,:,:,:);
    enc_stim  = stimval(opval==1);
    ret_phase = phase(opval==0,:,:,:);
    ret_stim  = stimval(opval==0);
    
    % cycle through every trial
    for i = 1 : size(enc_phase,1)
        for j = 1 : size(ret_phase,1)

            % if stimulus is not recalled, skip
            if ret_stim(j) == 0; continue; end
            
            % extract signal at encoding
            X = squeeze(enc_phase(i,:,:,freq.time>0&freq.time<=1));
                 
            % extract signal at retrieval
            Y = squeeze(ret_phase(j,:,:,freq.time>=0.5&freq.time<=1.5));
            
%             % cycle through each time window of retrieval
%             for t = 1 : size(ret_phase,4)-size(X,2)
%                         
%                 % get retieval signal for this trial
%                 Y(:,:,t) = squeeze(ret_phase(j,:,:,t:t+100));
%             end
            
            % get single-trial phase-locking 
            obs_model(count,:) = abs(mean(exp(1i.*(X-Y)),2));   
                
            % get prediction
            hyp_model(count,1) = enc_stim(i)==ret_stim(j);
            count = count + 1;        
        end
        if mod(i,10) == 0
            fprintf('Computing phase similarity... %d%% complete...\n',round((i/size(enc_phase,1))*100));
        end
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
