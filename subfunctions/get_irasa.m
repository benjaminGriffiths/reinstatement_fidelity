function [freq] = get_irasa(data,operation,modality)

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
freq.time           = 1;
freq.dimord         = 'rpt_chan_freq_time';
freq.trialinfo      = data.trialinfo;
freq.cfg            = [];

% restrict data to pre- and post-power
predata     = ft_selectdata(struct('latency',[-1 0]),data);
postdata    = ft_selectdata(struct('latency',[0.5 1.5]),data);

% extract variables to avoid overhead comm.
sig  = [predata.trial postdata.trial];
fs   = 100;
fr   = [5 25];

% predefine osci and frac results
osci = zeros(numel(sig),numel(data.label),26); 
fhz  = zeros(numel(sig),26);
frac = osci;

% cycle through every trial
fprintf('calculating irasa...\n')
parfor trl = 1 : numel(sig)
        
    % run irasa
    spec = amri_sig_fractal(sig{trl}',fs,...
                            'frange',fr);

    % record results
    osci(trl,:,:) = spec.osci';
    frac(trl,:,:) = spec.frac';
    fhz(trl,:)    = spec.freq;
end

% define derived freqs
foi = fhz(1,:);

% prepare to extract slop
fprintf('estimating intercept and slope...\n')
offset = zeros(size(frac,1)*size(frac,2),1);
slope = offset;
aucf  = offset;
auco  = offset;

% collapse power
frac = permute(frac,[3 1 2]);
frac = frac(:,:);
osc2 = permute(osci,[3 1 2]);
osc2 = osc2(:,:);

% get log-transformed data
X = log10(foi);
Y = log10(frac);   

% cycle through every point
parfor i = 1 : size(Y,2)
    
    % run linear model
    B = fitlm(X,squeeze(Y(:,i)));

    % store offset and slope
    offset(i,1) = B.Coefficients.tStat(1);
    slope(i,1) = B.Coefficients.tStat(2); 
    aucf(i,1) = trapz(foi,frac(:,i));
    auco(i,1) = trapz(foi,osc2(:,i));
end

% reshape output
offset = reshape(offset,[size(frac,1) size(frac,2)]);
slope  = reshape(slope,[size(frac,1) size(frac,2)]);

% add in freq
fprintf('packaging data...\n')
freq.powspctrm  = osci(numel(data.trial)+1:end,:,:);
freq.powslope   = slope(numel(data.trial)+1:end,:);
freq.powoffset  = offset(numel(data.trial)+1:end,:);
freq.powfracauc = aucf(numel(data.trial)+1:end,:);
freq.powosciauc = auco(numel(data.trial)+1:end,:);
freq.prepow     = osci(1:numel(data.trial),:,:);
freq.preslope   = slope(1:numel(data.trial),:);
freq.preoffset  = offset(1:numel(data.trial),:);
freq.preauc     = aucf(1:numel(data.trial),:);
freq.preosciauc = auco(1:numel(data.trial),:);
freq.freq       = foi;
fprintf('done...\n')
