function freq = get_phasepowercorr(data,operation,modality)

% define time/frequency resolution
if strcmpi(operation,'encoding')
    toi = data.time{1}(data.time{1}>=-1&data.time{1}<=1.5);
    foi = 8:30;
else    
    toi = data.time{1}(data.time{1}>=0&data.time{1}<=1.5);
    foi = 8:30;
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

% find nearest frequency to 8Hz
idx = knnsearch(freq.freq',8);

% extract phase and power
phase = angle(freq.fourierspctrm(:,:,idx,:));
power = abs(freq.fourierspctrm);

% ditch large data from freq structure
freq            = rmfield(freq,{'fourierspctrm','cumtapcnt','cfg'});
freq.dimord     = 'chan_freq_time';

% extract stim and operation values
stimval = nan(numel(freq.trialinfo),1);
opval = nan(numel(freq.trialinfo),1);
for trl = 1 : numel(freq.trialinfo)
    stimval(trl,1) = freq.trialinfo{trl}.stimulus_value;
    opval(trl,1) = strcmpi(freq.trialinfo{trl}.operation,'encoding');
end

% if is retrieval, limit power spectrum to save ram
if strcmpi(operation,'retrieval'); power = power(opval==0,:,:,freq.time>=0.5&freq.time<=1.5); end
    
% get time-averaged mean and standard deviation of power for each channel and frequency
avg_pow = mean(mean(power,4),1);
std_pow = std(mean(power,4),[],1);

% replicate matrices to match freq.powspctrm
avg_pow = repmat(avg_pow,[size(power,1) 1 1 size(power,4)]);
std_pow = repmat(std_pow,[size(power,1) 1 1 size(power,4)]);

% z-transform power
power = (power - avg_pow) ./ std_pow;
clear avg_pow std_pow

% change method based on whether encoding-encoding similarity or
% encoding-retreival similarity
if strcmpi(operation,'encoding')
       
    % get baseline difference in power
    baseidx = freq.time>=-1&freq.time<=0;
    postidx = freq.time>=0.5&freq.time<=1.5;
    power = nanmean(nanmean(power(:,:,:,postidx),4) - nanmean(power(:,:,:,baseidx),4),3);
    
    % predefine vectors for power and phase data
    phase_sim = nan(size(phase,1),size(phase,2));
    power_val = nan(size(phase,1),size(phase,2));
    
    % cycle through every trial
    for i = 1 : size(phase,1)
                
        % predefine observed and predicted model
        obs_model   = nan(size(phase,1),size(phase,2));
        hyp_model   = nan(size(phase,1),1);
        
        for j = 1 : size(phase,1)

            % skip matching trial
            if i == j; continue; end
            
            % extract signal
            X = squeeze(phase(i,:,:,:));
            Y = squeeze(phase(j,:,:,:));

            % get single-trial phase-locking 
            obs_model(j,:) = mean(cos(X).*cos(Y) + sin(X).*sin(Y),2);   

            % get prediction
            hyp_model(j,1) = stimval(i)==stimval(j);       
        end
            
        % remove nan rows
        hyp_model(all(isnan(obs_model),2),:) = [];
        obs_model(all(isnan(obs_model),2),:) = [];
        
        % cycle through each channel
        for chan = 1 : size(obs_model,2)
            
            % get similarity value for trial
            phase_sim(i,chan) = atanh(corr(obs_model(:,chan),hyp_model));
        end
        
        % get power value for trial
        power_val(i,:) = power(i,:);
        
        if mod(i,10) == 0
            fprintf('Computing phase similarity... %d%% complete...\n',round((i/size(phase,1))*100));
        end
    end
    
else
    
    % split encoding and retrieval data
    enc_phase = phase(opval==1,:,:,freq.time>0&freq.time<=1);
    enc_stim  = stimval(opval==1);
    ret_phase = phase(opval==0,:,:,freq.time>0.5&freq.time<=1.5);
    ret_stim  = stimval(opval==0);
    
    % get post-stim power
    ret_power = nanmean(nanmean(power,4),3);
    
    % predefine vectors for power and phase data
    phase_sim = nan(size(ret_phase,1),size(ret_phase,2));
    power_val = nan(size(ret_phase,1),size(ret_phase,2));
    
    % cycle through every trial
    for i = 1 : size(ret_phase,1)
                        
        % if stimulus is not recalled, skip
        if ret_stim(i) == 0; continue; end

        % predefine observed and predicted model
        obs_model   = nan(size(enc_phase,1),size(enc_phase,2));
        hyp_model   = nan(size(enc_phase,1),1);
        count       = 1;
        
        for j = 1 : size(enc_phase,1)
            
            % extract signal at encoding
            X = squeeze(enc_phase(j,:,:,:));
                 
            % extract signal at retrieval
            Y = squeeze(ret_phase(i,:,:,:));
            
            % get single-trial phase-locking 
            obs_model(count,:) = abs(mean(exp(1i.*(X-Y)),2));   
                
            % get prediction
            hyp_model(count,1) = enc_stim(j)==ret_stim(i);
            count = count + 1;        
        end
                
        % cycle through each channel
        for chan = 1 : size(obs_model,2)
            
            % get similarity value for trial
            phase_sim(i,chan) = atanh(corr(obs_model(:,chan),hyp_model));
        end
        
        % get power value for trial
        power_val(i,:) = ret_power(i,:);
        
        if mod(i,10) == 0
            fprintf('Computing phase similarity... %d%% complete...\n',round((i/size(ret_phase,1))*100));
        end
    end
    
    % remove nan columns
    power_val(all(isnan(power_val),2),:) = [];
    phase_sim(all(isnan(phase_sim),2),:) = [];
end
    
% predefine correlation vector
freq.powspctrm = nan(size(power_val,2),1);
freq.time       = 1;
freq.freq       = 1;

% cycle through each channel
for chan = 1 : size(power_val,2)
    
    % correlate observed and predicted model
    freq.powspctrm(chan,1) = atanh(corr(power_val(:,chan),phase_sim(:,chan)));
end
