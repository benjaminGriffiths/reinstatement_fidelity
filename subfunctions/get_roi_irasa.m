function [freq] = get_roi_irasa(data,roi,operation,modality)

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
freq.label          = {{'dummy'}};
freq.time           = 1;
freq.dimord         = 'rpt_chan_freq_time';
freq.cfg            = [];

% restrict data to pre- and post-power
A = ft_selectdata(struct('latency',[-1 0],'channels',{roi.label}),data);
B = ft_selectdata(struct('latency',[0.5 1.5],'channels',{roi.label}),data);

% extrat signal
Asig = reshape(cell2mat(A.trial),[size(A.trial{1},1) size(A.trial{1},2) numel(A.trial)]);
Bsig = reshape(cell2mat(B.trial),[size(B.trial{1},1) size(B.trial{1},2) numel(B.trial)]);

% extract variables to avoid overhead comm.
fs   = 100;
fr   = [5 25];

% predefine osci and frac results
osci = zeros(size(Bsig,3),26,2); 
fhz  = zeros(size(Bsig,3),26);
frac = osci;

% cycle through every trial
fprintf('calculating irasa...\n');
parfor trl = 1 : size(Bsig,3)
        
    % run irasa on post
    Bspec = amri_sig_avgfractal(Bsig(:,:,trl)',fs,...
                            'frange',fr);

    % run irasa on pre
    Aspec = amri_sig_avgfractal(Asig(:,:,trl)',fs,...
                            'frange',fr);
        
    % record results
    osci(trl,:,:) = [Bspec.osci Aspec.osci];
    frac(trl,:,:) = [Bspec.frac Aspec.frac];
    fhz(trl,:)  = Bspec.freq;    
end

% define derived freqs
foi = fhz(1,:);

% recode trialinfo
data = recode_trlinfo(data);

% add in freq
fprintf('packaging data...\n')
freq.powspctrm  = permute(osci,[1 3 2]);
freq.fractal    = permute(frac,[1 3 2]);
freq.freq       = foi;
freq.trialinfo  = data.trialinfo;
fprintf('done...\n')
