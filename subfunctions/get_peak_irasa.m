function [freq] = get_peak_irasa(data)

% check variables
fprintf('\npreparing data for irasa...\n')

% create freq structure
freq                = [];
freq.label          = data.label;
freq.time           = 1;
freq.dimord         = 'chan_freq_time';
freq.cfg            = [];

% extract signal
sig = reshape(cell2mat(data.trial),[size(data.trial{1},1) size(data.trial{1},2) numel(data.trial)]);
sig = permute(sig,[3 2 1]);

% extract variables to avoid overhead comm.
fs   = 100;
fr   = [5 25];

% predefine osci and frac results
osci = zeros(size(sig,3),205); 
fhz  = zeros(size(sig,3),205);
frac = osci;

% cycle through every trial
fprintf('calculating irasa...\n');
parfor chan = 1 : size(sig,3)
        
    % run irasa on post
    spec = amri_sig_avgfractal(sig(:,:,chan)',fs,...
                            'frange',fr);
                      
    % record results
    osci(chan,:) = spec.osci;
    frac(chan,:) = spec.frac;
    fhz(chan,:)  = spec.freq;    
end

% define derived freqs
foi = fhz(1,:);

% add in freq
fprintf('packaging data...\n')
freq.powspctrm  = osci;
freq.fractal    = frac;
freq.freq       = foi;
fprintf('done...\n')
