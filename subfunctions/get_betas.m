function [betas,h] = get_betas(freq)

% recode trialinfo
fprintf('preparing data for multi-regression...\n')
freq = recode_trlinfo(freq);

% extract parameters
pow     = freq.powspctrm;
slope   = freq.slope;

% find nans in power
nandx = any(isnan(pow))';

% create ERD outcome regressor
Y = double(1:size(pow,1) <= size(pow,1)/2)';

% create predictor matrix
X           = ones(size(Y,1),size(pow,2),3);
X(:,:,2)    = zscore(pow);
X(:,:,3)    = zscore(slope);

% predefine matrices for loop
b = zeros(2,size(X,2),size(X,3)); % type x channel x predictor
r = zeros(2,size(X,2));

% ---- ERD Regression ---- %
% cycle through encoding and retrieval
fprintf('running ERD multi-regression...\n')
for i = 1 : 2
    
    % restrict trials to task
    Yi = Y(freq.trialinfo(:,5) == i-1,1);
    Xi = X(freq.trialinfo(:,5) == i-1,:,:);

    % cycle through each channel
    for chan = 1 : size(pow,2)
        
        % check if nan
        if nandx(chan) == 1
            continue;
        end
        
        % restrict channels
        Yj = Yi;
        Xj = squeeze(Xi(:,chan,:));
        
        % test colinearity amongst predictors
        r(i,chan) = corr(Xj(:,2),Xj(:,3));
        
        % add constant
        b(i,chan,:) = Xj\Yj;
    end
end

% package result
fprintf('packaging ERD result...\n')
freq.dimord         = 'chan_freq_time';
betas               = repmat({rmfield(freq,{'powspctrm','slope','trialinfo'})},[4 1]);
betas{1}.powspctrm  = squeeze(b(2,:,2)); % encoding ERD power
betas{1}.slope      = squeeze(b(2,:,3)); % encoding ERD slope
betas{2}.powspctrm  = squeeze(b(1,:,2)); % retrieval ERD power
betas{2}.slope      = squeeze(b(1,:,3)); % retrieval ERD slope

% prepare plot
fprintf('plotting...\n') % design matrx, regressors, correlation coefficent
h = figure('position',[100 100 700 700]); hold on

% plot encoding
Xi = X(freq.trialinfo(:,5) == 1,:,:);
subplot(3,3,1); imagesc(zscore(squeeze(nanmean(Xi,2)))); ylabel('encoding')
title('design matrix'); set(gca,'xtick',1:3,'xticklabel',{'offset','power','slope'},'ytick',[]);
subplot(3,3,2); boxplot(squeeze(b(2,:,2:3))); yline(0);
title('resulting regressors'); set(gca,'xtick',1:2,'xticklabel',{'power','slope'});
subplot(3,3,3); histogram(r(2,:)); 
title('regressor collineratiry'); xlim([-1 1])

% plot retrieval
Xi = X(freq.trialinfo(:,5) == 0,:,:);
subplot(3,3,4); imagesc(zscore(squeeze(nanmean(Xi,2)))); ylabel('retrieval')
title('design matrix'); set(gca,'xtick',1:3,'xticklabel',{'offset','power','slope'},'ytick',[]);
subplot(3,3,5); boxplot(squeeze(b(1,:,2:3))); yline(0);
title('resulting regressors'); set(gca,'xtick',1:2,'xticklabel',{'power','slope'});
subplot(3,3,6); histogram(r(1,:)); 
title('regressor collineratiry'); xlim([-1 1])

% tidy workspace
clear Xi Yi Xj Yj b r i chan

% ---- RSE Regression ---- %
% extract trialinfo
trlinfo = freq.trialinfo;

% create memory model (cutting out prepow outcomes)
fprintf('running RSE multi-regression...\n')
Y = trlinfo(:,1);

% restrict trials to retrieval
Yi = Y(trlinfo(:,5) == 0,1);
Xi = X(trlinfo(:,5) == 0,:,:);

% drop pre-stim outcomes
Yi = Yi(1:size(Yi,1)/2,1);

% reshape X to include pre- and post- power as seperate regressors
Xi = cat(3,Xi(1:size(Yi,1),:,:),Xi(size(Yi,1)+1:end,:,2:3));

% predefine matrices for loop
b = zeros(size(Xi,2),size(Xi,3));
r = zeros(size(Xi,2),5,5);

% cycle through each channel
for chan = 1 : size(pow,2)
        
    % check if nan
    if nandx(chan) == 1
        continue;
    end

    % restrict channels
    Yj = Yi;
    Xj = squeeze(Xi(:,chan,:));

    % drop trials with no response
    Xj(isnan(Yj),:) = [];
    Yj(isnan(Yj),:) = [];
    
    % test colinearity amongst predictors
    r(chan,:,:) = corrcoef(Xj);

    % add constant
    b(chan,:) = Xj\Yj;
end

% package result
fprintf('packaging RSE result...\n')
betas{3}.powspctrm  = squeeze(b(:,2)); % poststim RSE power
betas{3}.slope      = squeeze(b(:,3)); % poststim RSE slope
betas{4}.powspctrm  = squeeze(b(:,4)); % prestim RSE power
betas{4}.slope      = squeeze(b(:,5)); % prestim RSE slope

% plot 
fprintf('plotting...\n') % plot encoding
subplot(3,3,7); imagesc(zscore(squeeze(nanmean(Xi,2)))); ylabel('RSE')
title('design matrix'); set(gca,'xtick',1:5,'xticklabel',{'offset','postpow','postslope','prepow','preslope'},'ytick',[]);
subplot(3,3,8); boxplot(squeeze(b(:,2:5))); yline(0);
title('resulting regressors'); set(gca,'xtick',1:4,'xticklabel',{'postpow','postslope','prepow','preslope'});
subplot(3,3,9); imagesc(squeeze(nanmean(r(:,2:5,2:5),1))); 
title('regressor collineratiry'); caxis([-1 1]);
fprintf('done...\n')

