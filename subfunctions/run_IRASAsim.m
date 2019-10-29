function run_IRASAsim

% --- INPUTS ---- %
operation = 'encoding';
modality = 'visual';
mask     = 1;
% -------------------%

% define root directory
if ispc;    dir_git  = 'E:/bjg335/projects/reinstatement_fidelity/';
            dir_bids = 'Y:/projects/reinstatement_fidelity/';
            dir_tool = 'Y:/projects/general/';
else;       dir_bids = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/bids_data/';
            dir_tool = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/general/';
end

% add developer rsa toolbox
addpath([dir_tool,'rsatoolbox-develop'])
addpath([dir_tool,'wen_spectral'])
n_subj      = 21;                                           % number of subjects

% load in similarity index
load([dir_bids,'derivatives/group/rsa-correlation/group_task-all_fmri-si'],'rsa_vec')
ssi = rsa_vec; % rename

% add in audio data
load([dir_bids,'derivatives/group/rsa-correlation-audio/group_task-all_fmri-si'],'rsa_vec')
ssi = cat(2,ssi,permute(rsa_vec,[1 3 2])); clear rsa_vec

% load mean bold
load([dir_bids,'derivatives/group/rsa-correlation/group_task-all_fmri-meanbold'],'mean_bold')
bold = mean_bold; % rename

% add in audio data
load([dir_bids,'derivatives/group/rsa-correlation-audio/group_task-all_fmri-meanbold'],'mean_bold')
bold = cat(2,bold,permute(mean_bold,[1 3 2])); clear mean_bold
   
% predefine structure for correlation values
grand_freq = repmat({struct('dimord','subj_chan_freq_time',...
                    'freq',8,...
                    'label',{{'dummy'}},...
                    'time',1,...
                    'cfg',[])},[3 1]); tic
                
% load mask details
load([dir_bids,'derivatives/group/rsa-correlation/masks.mat'],'mask_roi')

% cycle through each subject
for subj = 1 : n_subj

    % define subject data directory
    dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];    

    % load in raw data
    fprintf('\nloading sub-%02.0f data...\n',subj);
    load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'],'source')  
    
    % predefine conditional arrays to include all trials
    operation_to_include = zeros(numel(source.trial),1);
    modality_to_include  = zeros(numel(source.trial),1);

    % cycle through each trial
    for trl = 1 : numel(source.trial)
        operation_to_include(trl) = strcmpi(source.trialinfo{trl}.operation,operation);
        modality_to_include(trl) = strcmpi(source.trialinfo{trl}.modality,modality);
    end

    % select data
    cfg         = [];
    cfg.trials  = operation_to_include == 1 & modality_to_include == 1;
    cfg.channel = mask_roi{1}.label;
    cfg.latency = [0.5 1.5];
    source      = ft_selectdata(cfg,source);
    
    % get trial numbers
    source      = recode_trlinfo(source);
    trl_nums    = source.trialinfo(:,3);
    
    % extract bold for these trials
    bold_tmp = squeeze(bold(subj,mask,trl_nums));

    % fix numbers
    trl_nums(trl_nums>48 & trl_nums<=96)    = trl_nums(trl_nums>48 & trl_nums<=96) - 48;
    trl_nums(trl_nums>96 & trl_nums<=144)   = trl_nums(trl_nums>96 & trl_nums<=144) - 48;
    trl_nums(trl_nums>144 & trl_nums<=192)  = trl_nums(trl_nums>144 & trl_nums<=192) - 96;       

    % extract ssi
    ssi_tmp = squeeze(ssi(subj,mask,trl_nums));
    vals    = quantile(ssi_tmp,3);
    ssi_bin = nan(size(ssi_tmp));
    ssi_bin(ssi_tmp<=vals(1)) = 1;
    ssi_bin(ssi_tmp<=vals(2) & ssi_tmp>vals(1)) = 2;
    ssi_bin(ssi_tmp<=vals(3) & ssi_tmp>vals(2)) = 3;
    ssi_bin(ssi_tmp>vals(3)) = 4;
    
    % extract data
    dat = reshape(cell2mat(source.trial),[size(source.trial{1},1) size(source.trial{1},2) numel(source.trial)]);
    dat = permute(dat,[2 1 3]);
    
    % get binned bold
    bold_bin    = zeros(4,1);
    bold_bin(1) = mean(bold_tmp(ssi_bin==1));
    bold_bin(2) = mean(bold_tmp(ssi_bin==2));
    bold_bin(3) = mean(bold_tmp(ssi_bin==3));
    bold_bin(4) = mean(bold_tmp(ssi_bin==4));
    
    % get binned confidence
    bold_con    = zeros(4,1);
    bold_con(1) = nanmean(source.trialinfo(ssi_bin==1,2));
    bold_con(2) = nanmean(source.trialinfo(ssi_bin==2,2));
    bold_con(3) = nanmean(source.trialinfo(ssi_bin==3,2));
    bold_con(4) = nanmean(source.trialinfo(ssi_bin==4,2));
    
    % predefine output vecs
    pow = zeros(4,1);
    ofs = zeros(4,1);
    slp = zeros(4,1);
    
    % cycle through each bin
    for i = 1 : 4
        
        % get irasa of binned psd
        tmp = dat(:,:,ssi_bin==i);
        output = amri_sig_avgfractal(tmp(:,:),100,'frange',[5 25]);
    
        % get alpha/beta power
        pow(i,1) = mean(output.osci(output.freq>=8));
        
        % log transform fractal
        logF = log10(output.freq');
        logP = squeeze(log10(output.frac))';

        % get linfit
        B = fitlm(logF,logP); %post

        % get slope and intercept difference
        ofs(i,1) = B.Coefficients.tStat(1);
        slp(i,1) = B.Coefficients.tStat(2);
    end

    % create regressor matrix
    X      = zeros(4,3);
    X(:,1) = [1 2 3 4];
    X(:,2) = zscore(pow);
    X(:,3) = zscore(slp);
    
    % create table
    tbl = array2table(X,'VariableNames',{'rsa','pow','slp'});

    % define predictors
    preds = {'pow','slp'};

    % run model
    B = fitlm(tbl,'ResponseVar','rsa',...
                  'PredictorVars',preds);
      
    % store outputs
    regweights(subj,:) = B.Coefficients.tStat(2:end);
    
    % update user
    fprintf('sub-%02.0f complete...\n',subj)
end

% save data
save([dir_bids,'derivatives/group/rsa-correlation/group_task-all_comb-freq.mat'],'grand_freq'); 
