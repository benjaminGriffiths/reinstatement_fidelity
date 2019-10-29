function freq = get_roi_fooof(data,roi,modality,operation)

% check variables
fprintf('organising data...\n')

% select channels of interest
% get modalities and tasks
ismodality    = zeros(numel(data.trialinfo),1);
isoperation  = zeros(numel(data.trialinfo),1);

% cycle through trials
for trl = 1 : numel(data.trialinfo)
    ismodality(trl) = strcmpi(data.trialinfo{trl}.modality,modality);
    isoperation(trl) = strcmpi(data.trialinfo{trl}.operation,operation);
end

% recode trialinfo
data = recode_trlinfo(data);

% predefine freq structure
freq = struct('time',1,...
              'freq',8,...
              'label',{{'dummy'}},...
              'dimord','chan_freq_time');

% define trials  and channels of interest
trls  	= ismodality & isoperation;
chans   = ismember(data.label,roi.label);
times   = data.time{1} >= 0.5 & data.time{1} <= 1.5;
prestim = data.time{1} >= -1 & data.time{1} <= 0;

% define key variables
ntrl    = numel(data.trial(trls));
nchan   = numel(data.label);
ntime   = numel(data.time{1});

% update user
fprintf('calculating PSD...\n');

% extract post-stim signal
signal = permute(reshape(cell2mat(data.trial(trls)),[nchan,ntime,ntrl]),[2 1 3]);
signal = signal(times,chans,:);

% get psd
[A,freqs] = pwelch(signal(:,:),50,25,[],100);

% extract pre-stim signal
signal = permute(reshape(cell2mat(data.trial(trls)),[nchan,ntime,ntrl]),[2 1 3]);
signal = signal(prestim,chans,:);

% get psd
B = pwelch(signal(:,:),50,25,[],100);

% reshape psd and average over channels
A = squeeze(mean(reshape(A,[size(A,1),size(signal,2),size(signal,3)]),2));
B = squeeze(mean(reshape(B,[size(B,1),size(signal,2),size(signal,3)]),2));

% fit fooof
[Ac,As,Ap] = get_fooof_diff(A,freqs);
[Bc,Bs,Bp] = get_fooof_diff(B,freqs);

% package post-stim data
freq.powspctrm = Ac(:,2);
freq.bckgrnd   = Ac(:,1);
freq.intercept = As(:,1);
freq.slope     = As(:,2);
freq.peakfreq  = Ap(:,1);
freq.peakpow   = Ap(:,2);
freq.trialinfo = data.trialinfo(trls,:);

% package post-stim data
freq.prespctrm = Bc(:,2);
freq.prebck    = Bc(:,1);
freq.preintcpt = Bs(:,1);
freq.preslope  = Bs(:,2);
freq.prepkfrq  = Bp(:,1);
freq.prepkpow  = Bp(:,2);

% package difference data
freq.erdspctrm = Ac(:,2) - Bc(:,2);
freq.erdslope  = Ac(:,1) - Bc(:,1);
freq.erdintcpt = As(:,1) - Bs(:,1);
freq.erdslope  = As(:,2) - Bs(:,2);
freq.erdpkfrq  = Ap(:,1) - Bp(:,1);
freq.erdpkpow  = Ap(:,2) - Bp(:,2);

end

function [auc,slp,pks] = get_fooof_diff(A,freqs)

% predefine output data
fprintf('Calculating FOOOF...\n')
auc = nan(size(A,2),2);
slp = nan(size(A,2),2);
pks = nan(size(A,2),2);

% cycle through each trl
parfor i = 1 : size(A,2)

    % initialize FOOOF object
    fm = py.fooof.FOOOF([1 8],...    % peak width
                        8,...        % n peaks
                        0,...        % min amp.
                        2,...        % peak thr.
                        'fixed',...  % knee
                        false);      % verbose

    % convert inputs
    pyF   = py.numpy.array(freqs');
    f_range = py.list([5 30]);
    
    % extract power
    pyA = py.numpy.array(A(:,i)');
        
    % run FOOOF fit on A
    fm.fit(pyF, pyA, f_range)

    % extract outputs
    fitA = fm.get_results();
    fitA = fooof_unpack_results(fitA);
    modA = fooof_get_model(fm);
    for field = fieldnames(modA)'
        fitA.(field{1}) = modA.(field{1});
    end

    % get difference in oscillatory spectrum
    fidx = (fitA.freqs>=7 & fitA.freqs<=14);
    oscA = fitA.fooofed_spectrum(fidx) - fitA.bg_fit(fidx);    
    
    % get area under slope
    auc(i,:) = [trapz(fitA.bg_fit),trapz(oscA)];    
    
    % get peak alpha
    pp = fitA.peak_params(fitA.peak_params(:,1) <= 14,1:2);
    
    % if alpha peak exists
    if ~isempty(pp)
        [~,j] = max(pp(:,1));
        pks(i,:) = pp(j,:);
    end
    
    % get slope function
    slp(i,:) = fitA.background_params;
end
end
