function [freq,signal,h] = get_condition_fooof(data)

% check variables
fprintf('organising data...\n')

% get modalities and tasks
isvisual    = zeros(numel(data.trialinfo),1);
isencoding  = zeros(numel(data.trialinfo),1);

% cycle through trials
for trl = 1 : numel(data.trialinfo)
    isvisual(trl) = strcmpi(data.trialinfo{trl}.modality,'visual');
    isencoding(trl) = strcmpi(data.trialinfo{trl}.operation,'encoding');
end

% recode trialinfo
data    = recode_trlinfo(data);
trial   = data.trial;

% % filter
% fprintf('filtering data...\n')
% parfor trl = 1 : numel(data.trial)
%     dat = trial{trl};
% 	dat = ft_preproc_bandstopfilter(dat,100,[15 17],4,'but','twopass','no');
% 	trial{trl} = ft_preproc_bandstopfilter(dat,100,[30 34],4,'but','twopass','no');
% end
% data.trial = trial;

% predefine freq structure
freq = repmat({struct('time',1,...
              'freq',8,...
              'label',{data.label},...
              'powspctrm',[],...
              'slope',[],...
              'dimord','chan_freq_time')},[2 3]);

% predefine signal output
signal = cell(2,3);
          
% cycle through conditions
for i = 1 : 2
    for j = 1 : 2

        % define trials of interest
        trls = isvisual == (i==1) & isencoding == (j==1);
        
        % calculate erd
        [A,B,freqs] = get_erd_psd(data,trls);

        % fit fooof
        [auc,signal{i,j},params] = get_fooof_diff(A,B,freqs);

        % package data
        freq{i,j}.powspctrm = auc(:,2);
        freq{i,j}.slope     = auc(:,1);
        freq{i,j}.offset    = params(:,1);
        freq{i,j}.beta      = params(:,2);
        
        % if retrieval
        if j == 2
            
            % calculate rse psd
            [A,B,freqs] = get_rse_psd(data,trls);

            % fit fooof
            [auc,signal{i,j+1},params] = get_fooof_diff(A,B,freqs);

            % package data
            freq{i,j+1}.powspctrm   = auc(:,2);
            freq{i,j+1}.slope       = auc(:,1);
            freq{i,j+1}.offset      = params(:,1);
            freq{i,j+1}.beta        = params(:,2);
        end
    end
end

% tidy workspace
freq = freq';
freq = freq(:);
fprintf('done...\n')

% plot if requested
if nargout > 1
    
    % create figure
    h = figure('visible','off'); hold on
    count = 1;
     
    % define label names
    cN = {'Enc. ERD','Ret. ERD','Ret. Succ.'};
    
    % cycle through conditions
    for i = 1 : 3
        
        % plot spectrum
        subplot(3,3,count); hold on; count = count + 1;
        sy = signal{1,i};
        fx = signal{1,i}(:,5,1);
        shadedErrorBar(fx,mean(sy(:,3,:),3),sem(sy(:,3,:),3),{'color',[0.7 0.5 0.4]})
        shadedErrorBar(fx,mean(sy(:,4,:),3),sem(sy(:,4,:),3),{'color',[0.4 0.5 0.7]})
        if i == 1; title('Combined'); end
        ylabel(cN{i})
        
        % plot slope
        subplot(3,3,count); hold on; count = count + 1;
        sy = signal{1,i};
        fx = signal{1,i}(:,5,1);
        shadedErrorBar(fx,mean(sy(:,1,:),3),sem(sy(:,1,:),3),{'color',[0.7 0.5 0.4]})
        shadedErrorBar(fx,mean(sy(:,2,:),3),sem(sy(:,2,:),3),{'color',[0.4 0.5 0.7]})
        if i == 1; title('Aperiodic'); end
        
        % plot oscillation
        subplot(3,3,count); hold on; count = count + 1;
        sy = signal{1,i};
        fx = signal{1,i}(:,5,1);
        shadedErrorBar(fx,mean(sy(:,3,:)-sy(:,1,:),3),sem(sy(:,3,:)-sy(:,1,:),3),{'color',[0.7 0.5 0.4]})
        shadedErrorBar(fx,mean(sy(:,4,:)-sy(:,2,:),3),sem(sy(:,4,:)-sy(:,2,:),3),{'color',[0.4 0.5 0.7]})
        if i == 1; title('Periodic'); end
    end    
end

end

function [A,B,freqs] = get_erd_psd(data,trls)

% update user
fprintf('calculating ERD FOOOF...\n');

% define key variables
ntrl    = numel(data.trial(trls));
nchan   = numel(data.label);
ntime   = numel(data.time{1});

% extract signal
signal = permute(reshape(cell2mat(data.trial(trls)),[nchan,ntime,ntrl]),[2 1 3]);

% get post-stim psd
postsig = signal(data.time{1}>=0.5&data.time{1}<=1.5,:,:);
[A,freqs] = pwelch(postsig(:,:),50,25,[],100);

% get pre-stim psd
presig  = signal(data.time{1}>=-1&data.time{1}<=0,:,:);
[B,~] = pwelch(presig(:,:),50,25,[],100);

% reshape psd and average over channels
A = mean(reshape(A,[size(A,1),size(signal,2),size(signal,3)]),3);
B = mean(reshape(B,[size(B,1),size(signal,2),size(signal,3)]),3);

end


function [A,B,freqs] = get_rse_psd(data,trls)

% update user
fprintf('calculating RSE PSD...\n');

% define key variables
ntrl    = numel(data.trial(trls));
nchan   = numel(data.label);
ntime   = numel(data.time{1});
memory  = data.trialinfo(trls,1);

% extract signal
signal = permute(reshape(cell2mat(data.trial(trls)),[nchan,ntime,ntrl]),[2 1 3]);

% get post-stim psd
hits        = signal(data.time{1}>=0.5&data.time{1}<=1.5,:,memory==1);
[A,freqs]   = pwelch(hits(:,:),50,25,[],100);

% get pre-stim psd
misses      = signal(data.time{1}>=0.5&data.time{1}<=1.5,:,memory==0);
[B,~]       = pwelch(misses(:,:),50,25,[],100);

% reshape psd and average over channels
A = mean(reshape(A,[size(A,1),size(signal,2),size(hits,3)]),3);
B = mean(reshape(B,[size(B,1),size(signal,2),size(misses,3)]),3);

end


function [auc,sig,params] = get_fooof_diff(A,B,freqs)

% predefine output data
auc = nan(size(A,2),2);
params = nan(size(A,2),2);
sig = nan(95,5,size(A,2));

% cycle through each trl
parfor i = 1 : size(A,2)

    % initialize FOOOF object
    fm = py.fooof.FOOOF([1 8],...    % peak width
                        inf,...        % n peaks
                        0,...        % min amp.
                        2,...        % peak thr.
                        'fixed',...  % knee
                        false);      % verbose

    % convert inputs
    pyF   = py.numpy.array(freqs');
    f_range = py.list([3 40]);
    
    % extract power
    pyA = py.numpy.array(A(:,i)');
    pyB = py.numpy.array(B(:,i)');
        
    % run FOOOF fit on A
    fm.fit(pyF, pyA, f_range)

    % extract outputs
    fitA = fm.get_results();
    fitA = fooof_unpack_results(fitA);
    modA = fooof_get_model(fm);
    for field = fieldnames(modA)'
        fitA.(field{1}) = modA.(field{1});
    end

    % run FOOOF fit on B
    fm.fit(pyF, pyB, f_range)

    % extract outputs
    fitB = fm.get_results();
    fitB = fooof_unpack_results(fitB);
    modB = fooof_get_model(fm);
    for field = fieldnames(modB)'
        fitB.(field{1}) = modB.(field{1});
    end
    
    % extract signal
    sig(:,:,i) = [fitA.bg_fit;fitB.bg_fit;fitA.fooofed_spectrum;fitB.fooofed_spectrum;fitB.freqs]';
    
    % get difference in aperiodic signal
    fracDiff = trapz(fitA.bg_fit) - trapz(fitB.bg_fit);
    
    % get difference in oscillatory spectrum
    fidx = (fitA.freqs>=6 & fitA.freqs<=14);
    oscA = fitA.fooofed_spectrum(fidx) - fitA.bg_fit(fidx);
    oscB = fitB.fooofed_spectrum(fidx) - fitB.bg_fit(fidx);
    
    % get area under oscillatory curve
    oscDiff = trapz(oscA) - trapz(oscB);
    
    % get area under slope
    auc(i,:) = [fracDiff,oscDiff];
    
    % extract params
    params(i,:) = fitA.background_params - fitB.background_params;
                   
end
end
