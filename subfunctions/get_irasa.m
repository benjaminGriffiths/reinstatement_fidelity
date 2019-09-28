function [freq] = get_irasa(data,operation,modality,contrast)

% check variables
fprintf('\npreparing data for irasa...\n')

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
config        = [];
config.trials = operation_to_include == 1 & modality_to_include == 1;
data          = ft_selectdata(config,data);

% create freq structure
freq                = [];
freq.label          = data.label;
freq.time           = [-0.5 1];
freq.dimord         = 'chan_freq_time';
freq.cfg            = [];

% if ERD contrast
if strcmpi(contrast,'erd')

    % restrict data to pre- and post-power
    A = ft_selectdata(struct('latency',[-1 0]),data);
    B = ft_selectdata(struct('latency',[0.5 1.5]),data);

    % extrat signal
    Asig = reshape(cell2mat(A.trial),[size(A.trial{1},1) size(A.trial{1},2) numel(A.trial)]);
    Asig = permute(Asig,[3 2 1]);
    Bsig = reshape(cell2mat(B.trial),[size(B.trial{1},1) size(B.trial{1},2) numel(B.trial)]);
    Bsig = permute(Bsig,[3 2 1]);
    
% if  RSE contrast    
else
    
    % recode trial info
    tmp = recode_trlinfo(data);
    
    % restrict data to pre- and post-power
    A = ft_selectdata(struct('latency',[0.5 1.5],'trials',tmp.trialinfo(:,1) == 0),data);
    B = ft_selectdata(struct('latency',[0.5 1.5],'trials',tmp.trialinfo(:,1) == 1),data);

    % extrat signal
    Asig = reshape(cell2mat(A.trial),[size(A.trial{1},1) size(A.trial{1},2) numel(A.trial)]);
    Asig = permute(Asig,[3 2 1]);
    Bsig = reshape(cell2mat(B.trial),[size(B.trial{1},1) size(B.trial{1},2) numel(B.trial)]);
    Bsig = permute(Bsig,[3 2 1]);
end

% extract variables to avoid overhead comm.
fs   = 100;
fr   = [3 25];

% predefine osci and frac results
osci = zeros(size(Bsig,2),29,2); 
fhz  = zeros(size(Bsig,2),29);
frac = osci;

% cycle through every trial
fprintf('calculating irasa...\n');
for chan = 1 : size(Bsig,3)
        
    % run irasa on post
    Bspec = amri_sig_avgfractal(Bsig(:,:,chan)',fs,...
                            'frange',fr);

    % run irasa on pre
    Aspec = amri_sig_avgfractal(Asig(:,:,chan)',fs,...
                            'frange',fr);
                        
    % record results
    osci(chan,:,:) = [Aspec.osci Bspec.osci];
    frac(chan,:,:) = [Aspec.frac Bspec.frac];
    fhz(chan,:)    = Aspec.freq;    
end

% define derived freqs
foi = fhz(1,:);

% add in freq
fprintf('packaging data...\n')
freq.powspctrm  = osci;
freq.fractal    = frac;
freq.freq       = foi;
fprintf('done...\n')
