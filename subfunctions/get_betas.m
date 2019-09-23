function betas = get_betas(freq,regtype)

% recode trialinfo
freq = recode_trlinfo(freq);

% extract parameters
pow     = freq.powspctrm;
slope   = freq.slope;

% switch based on regression type requested
switch regtype
    
    % get ERD design matrix
    case 'erd'        
        [X,Y] = prepERD(pow,slope);

    % get RSE design matrix
    case 'rse'
        [X,Y]   = prepRSE(pow,slope,freq.trialinfo);
end

% run regressor
b = runreg(X,Y);

% package result
betas = struct('cfg',[],...
               'freq',8,...
               'time',1,...
               'label',{freq.label},...
               'powspctrm',b(:,1),...
               'slope',b(:,2),...
               'dimord','chan_freq_time');

end

function b = runreg(X,Y)

% predefine matrices for loop
b = zeros(size(X,2),size(X,3)-1);

% cycle through each channel
for chan = 1 : size(X,2)

    % restrict channels
    Xj = squeeze(X(:,chan,:));

    % run regression
    bout        = Xj\Y;
    b(chan,:)   = bout(2:end);
end
end

function [X,Y] = prepERD(pow,slope)

% update user
fprintf('preparing ERD multi-regression...\n')

% create ERD outcome regressor
Y = double(1:size(pow,1) <= size(pow,1)/2)';

% create drift regressor
lindrf = zscore(1:size(Y,1)/2)';

% create predictor matrix
X           = ones(size(Y,1),size(pow,2),4)-0.5;
X(:,:,2)    = zscore(pow);
X(:,:,3)    = zscore(slope);
X(:,:,4)    = repmat(lindrf,[2 size(X,2)]);

end

function [X,Y]   = prepRSE(pow,slope,trlinfo)

% update user
fprintf('preparing RSE multi-regression...\n')

% get pre/post split
postidx = 1:size(trlinfo,1)/2;
preidx  = size(trlinfo,1)/2+1:size(trlinfo,1);

% create ERD outcome regressor
Y = trlinfo(postidx,1);

% extract confidence from trlinfo
conf = trlinfo(postidx,2);
conf = (conf - nanmean(conf)) ./ nanstd(conf);

% create predictor matrix
X           = ones(size(Y,1),size(pow,2),7)-0.5;
X(:,:,2)    = zscore(pow(postidx,:));
X(:,:,3)    = zscore(slope(postidx,:));
X(:,:,4)    = zscore(pow(preidx,:));
X(:,:,5)    = zscore(slope(preidx,:));
X(:,:,6)    = repmat(conf,[1 size(X,2)]);
X(:,:,7)    = repmat(zscore(1:size(Y,1))',[1 size(X,2)]);

% remove nans
nandx       = isnan(Y) | any(isnan(X(:,1,:)),3);
Y(nandx)    = [];
X(nandx,:,:) = [];

end
