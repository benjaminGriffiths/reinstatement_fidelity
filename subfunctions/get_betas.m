function [betas,h] = get_betas(freq,method)

% check inputs
fprintf('preparing data for multi-regression...\n')
if nargin == 1; method = 'linreg'; end

% recode trialinfo
freq = recode_trlinfo(freq);

% extract parameters
pow     = freq.powspctrm;
slope   = freq.slope;

% create ERD outcome regressor
Y = double(1:size(pow,1) <= size(pow,1)/2)';

% create predictor matrix
X           = ones(size(Y,1),size(pow,2),3);
X(:,:,1)    = zscore(pow);
X(:,:,2)    = zscore(slope);

% define linear regressors
lX = freq.trialinfo(:,[4 3]);

% predefine matrices for loop
b = zeros(2,size(X,2),size(X,3)); % type x channel x predictor
r = zeros(2,size(X,2));

% ---- ERD Regression ---- %
% cycle through encoding and retrieval
fprintf('running ERD multi-regression...\n')
for i = 1 : 2
    
    % add linear regressor
    X(:,:,3) = repmat(lX(:,i),[1 size(X,2)]);
    
    % restrict trials to task
    Yi = Y(freq.trialinfo(:,5) == i-1,1);
    Xi = X(freq.trialinfo(:,5) == i-1,:,:);

    % sort linear regressor
    Xs = tiedrank(Xi(:,1,3));
    Xi(:,:,3) = repmat(Xs,[1 size(Xi,2)]);
    
    % cycle through each channel
    for chan = 1 : size(pow,2)
        
        % restrict channels
        Yj = Yi;
        Xj = squeeze(Xi(:,chan,:));
        
        % test colinearity amongst predictors
        r(i,chan) = corr(Xj(:,1),Xj(:,2));
        
        % run regression
        switch method
            case 'linreg'
                B = fitlm(Xj,Yj);
                b(i,chan,:) = B.Coefficients.tStat(2:end);
            case 'logreg'
                [~,~,t]= mnrfit(Xj,Yj+1);
                b(i,chan,:) = -t.t(2:end);
            case 'robust'
                B = fitlm(Xj,Yj,'RobustOpts','on');
                b(i,chan,:) = B.Coefficients.tStat(2:end);
        end
    end
end

% package result
fprintf('packaging ERD result...\n')
freq.dimord         = 'chan_freq_time';
betas               = repmat({rmfield(freq,{'powspctrm','slope','trialinfo'})},[4 1]);
betas{1}.powspctrm  = squeeze(b(2,:,1))'; % encoding ERD power
betas{1}.slope      = squeeze(b(2,:,2))'; % encoding ERD slope
betas{2}.powspctrm  = squeeze(b(1,:,1))'; % retrieval ERD power
betas{2}.slope      = squeeze(b(1,:,2))'; % retrieval ERD slope

% prepare plot
fprintf('plotting...\n') % design matrx, regressors, correlation coefficent
h = figure('position',[100 100 700 700]); hold on

% plot encoding
Xi = X(freq.trialinfo(:,5) == 1,:,:);
subplot(3,3,1); imagesc(zscore(squeeze(nanmean(Xi,2)))); ylabel('encoding')
title('design matrix'); set(gca,'xtick',1:3,'xticklabel',{'power','slope','linear'},'ytick',[]);
subplot(3,3,2); boxplot(squeeze(b(2,:,:))); yline(0);
title('resulting regressors'); set(gca,'xtick',1:3,'xticklabel',{'power','slope','linear'});
subplot(3,3,3); histogram(r(2,:)); 
title('regressor collineratiry'); xlim([-1 1])

% plot retrieval
Xi = X(freq.trialinfo(:,5) == 0,:,:);
subplot(3,3,4); imagesc(zscore(squeeze(nanmean(Xi,2)))); ylabel('retrieval')
title('design matrix'); set(gca,'xtick',1:3,'xticklabel',{'power','slope','linear'},'ytick',[]);
subplot(3,3,5); boxplot(squeeze(b(1,:,:))); yline(0);
title('resulting regressors'); set(gca,'xtick',1:3,'xticklabel',{'power','slope','linear'});
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

% add linear regressor
X(:,:,3) = repmat(lX(:,1),[1 size(X,2)]);

% restrict trials to retrieval
Yi = Y(trlinfo(:,5) == 0,1);
Xi = X(trlinfo(:,5) == 0,:,:);

% sort linear regressor
Xs = tiedrank(Xi(:,1,3));
Xi(:,:,3) = repmat(Xs,[1 size(Xi,2)]);

% drop pre-stim outcomes
Yi = Yi(1:size(Yi,1)/2,1);

% reshape X to include pre- and post- power as seperate regressors
Xi = cat(3,Xi(1:size(Yi,1),:,1:2),Xi(size(Yi,1)+1:end,:,:));

% add confidence
conf = trlinfo(trlinfo(:,5) == 0,2);
Xi(:,:,6) = repmat(conf(1:size(Yi,1)),[1,size(Xi,2)]);

% identify trials with no response
nandx = any(isnan(cat(2,squeeze(Xi(:,1,:)),Yi)),2);

% predefine matrices for loop
b = zeros(size(Xi,2),size(Xi,3));
r = zeros(size(Xi,2),6,6);

% cycle through each channel
for chan = 1 : size(pow,2)
     
    % restrict channels
    Yj = Yi;
    Xj = squeeze(Xi(:,chan,:));

    % drop trials with no response
    Xj(nandx,:) = [];
    Yj(nandx,:) = [];
    
    % test colinearity amongst predictors
    r(chan,:,:) = corrcoef(Xj);

    % run regression
    switch method
        case 'linreg'
            B = fitlm(Xj,Yj);
            b(chan,:) = B.Coefficients.tStat(2:end);
        case 'logreg'
            [~,~,t]= mnrfit(Xj,Yj+1);
            b(chan,:) = -t.t(2:end);
        case 'robust'
            B = fitlm(Xj,Yj,'RobustOpts','on');
            b(chan,:) = B.Coefficients.tStat(2:end);
    end
end

% package result
fprintf('packaging RSE result...\n')
betas{3}.powspctrm  = squeeze(b(:,1)); % poststim RSE power
betas{3}.slope      = squeeze(b(:,2)); % poststim RSE slope
betas{4}.powspctrm  = squeeze(b(:,3)); % prestim RSE power
betas{4}.slope      = squeeze(b(:,4)); % prestim RSE slope

% plot 
fprintf('plotting...\n') % plot encoding
subplot(3,3,7); imagesc(zscore(squeeze(nanmean(Xi,2)))); ylabel('RSE')
title('design matrix'); set(gca,'xtick',1:6,'xticklabel',{'postpow','postslope','prepow','preslope','lin','conf'},'ytick',[]);
subplot(3,3,8); boxplot(b); yline(0);
title('resulting regressors'); set(gca,'xtick',1:6,'xticklabel',{'postpow','postslope','prepow','preslope','lin','conf'});
subplot(3,3,9); imagesc(squeeze(nanmean(r,1))); 
title('regressor collineratiry'); caxis([-1 1]);
fprintf('done...\n')

end
