function run_EEG_analysis(skipComputation)

% if skip argument not entered, set skip to true
if nargin == 0; skipComputation = true; end

% define root directory
if ispc;    dir_bids = 'Y:/projects/reinstatement_fidelity/bids_data/';
            dir_repos = 'E:/bjg335/projects/'; % repository directory
else;       dir_bids = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/reinstatement_fidelity/bids_data/';
            dir_repos = '/media/bjg335/rds-share-2018-hanslmas-memory/projects/git_clone/'; % repository directory
end

% add subfunctions
addpath([dir_repos,'reinstatement_fidelity/subfunctions'])
addpath(genpath([dir_repos,'fooof_mat/']))

% define number of subjects
n_subj = 21;
n_chan = 1817;

%% Get Power and Fractal for Each Participant
% get IRASA if requested
if ~skipComputation; tic

    % start parallel pool
    p = gcp('nocreate');
    if isempty(p)
        parpool(4);
    end
    
    % cycle through each subject
    for subj = 1 : n_subj

        % define subject data directory
        dir_subj = [dir_bids,'sourcedata/',sprintf('sub-%02.0f',subj),'/eeg/'];    

        % load in raw data
        fprintf('\nloading sub-%02.0f data...\n',subj);
        load([dir_subj,sprintf('sub-%02.0f',subj),'_task-rf_eeg-source.mat'],'source')  

        % get irasa
        freq{1} = get_irasa(source,'encoding','visual','erd');
        freq{2} = get_irasa(source,'retrieval','visual','erd');
        freq{3} = get_irasa(source,'retrieval','visual','rse');

        % save spectral outputs
        save(sprintf('%sderivatives/sub-%02.0f/eeg/sub-%02.0f_task-rf_eeg-irasa.mat',dir_bids,subj,subj),'freq')  
        fprintf('\nsub-%02.0f complete...\ntotal time elapsed: %1.1f minutes...\nestimated time remaining: %1.1f minutes...\n',subj,round(toc/60,1),round((ound(toc/60,1)/subj)*(n_subj-subj),1))

        % tidy up
        close all
        clear dir_subj subj freq spec h source
    end
    
    % close parallel pool
    delete(gcp);
end

%% Extract Slope and Peak Alpha Power
% predefine group freq structure
grand_freq = repmat({struct('cfg',[],'time',1,...
                    'freq',8,'label',{{}},...
                    'dimord','subj_chan_freq_time',...
                    'powspctrm',nan(n_subj,n_chan),...
                    'slope',nan(n_subj,n_chan))},[6 1]);

% cycle through each subject
for subj = 1 : n_subj
    
    % load in raw data
    load(sprintf('%sderivatives/sub-%02.0f/eeg/sub-%02.0f_task-rf_eeg-irasa.mat',dir_bids,subj,subj),'freq')  
        
    % plot freq spectrum for each output
    h = int_plotFreqSpec(freq);
    saveas(h,sprintf('%s/derivatives/group/figures/irasa/sub-%02.0f_irasa-visual',dir_bids,subj),'jpg')
    close all
        
    % cycle through each cell
    for i = 1 : numel(freq)
    
        % get alpha power difference
        alpha_idx = freq{i}.freq >= 7 & freq{i}.freq <= 25;
        grand_freq{i}.powspctrm(subj,:) = mean(freq{i}.powspctrm(:,alpha_idx,2) - freq{i}.powspctrm(:,alpha_idx,1),2);
        grand_freq{i}.label = freq{i}.label;
        
        % cycle through each channel
        for chan = 1 : size(freq{i}.fractal,1)
            
            % get fractal
            A = freq{i}.fractal(chan,:,1)';
            B = freq{i}.fractal(chan,:,2)';
            f = freq{i}.freq';
            
            % set intercept to zero
            A = detrend(A,0);
            B = detrend(B,0);
            
            % get slope
            grand_freq{i}.slope(subj,chan) = B\f - A\f;
        end
    end
    
    % update user    
    fprintf('sub-%02.0f done...\n',subj);
end

% load ROI data
load([dir_bids,'derivatives/group/eeg/roi.mat'],'roi')

% cycle through each data type
for i = 1 : numel(grand_freq)
    
    % convert datatype to source
    grand_freq{i,1} = int_convertData(grand_freq{i,1},roi);
end

% save data
save([dir_bids,'derivatives/group/eeg/group_task-all_pow-fooof.mat'],'grand_freq');

%% Run Statistics
% set seed
rng(1) 

% define names of analyses
metric = {'pow','slp'};
condition = {'vis_enc','vis_ret','vis_rse','aud_enc','aud_ret','aud_rse'};

% run statistics
[stat,tbl] = int_runStats(grand_freq,metric,[-1 0]);
disp(tbl)

% save data
save([dir_bids,'derivatives/group/eeg/group_task-all_pow-stat.mat'],'stat','tbl');

% extract subject activation in cluster
clus_tbl = int_extractSubjs(grand_freq,stat,tbl,metric,condition);

% write table
writetable(clus_tbl,[dir_bids,'derivatives/group/eeg/group_task-all_eeg-cluster.csv'],'Delimiter',',')
copyfile([dir_bids,'derivatives/group/eeg/group_task-all_eeg-cluster.csv'],[dir_repos,'reinstatement_fidelity/data/fig2_data/group_task-rf_eeg-cluster.csv'])

% cycle through each condition
for i = 1 : numel(condition)
    
    % get power brain map
    filename = [dir_bids,'derivatives/group/eeg/group_task-',condition{i},'_eeg-map.nii'];
    int_extractMap(stat{i},filename);
    
    % copy to git repo
    copyfile(filename,[dir_repos,'reinstatement_fidelity/data/fig2_data/group_task-',condition{i},'_eeg-map.nii'])    
end

end

% ----------------------------------------------------------------------- %
%                                                                         %
%    SUBFUNCTIONS                                                         %
%    ------------                                                         %
%       ------                                                            %
%                                                                         %
% ----------------------------------------------------------------------- %

function h = int_plotFreqSpec(freq)

% plot outputs
h = figure('visible','off','position',[100 100 800 400]); hold on; count = 1;
col = [0.7 0.4 0.3; 0.3 0.4 0.7; 0.4 0.5 0.2];
lab = {'ENC:','RET:','RSE:'};
for i = 1 : 3
    subplot(3,4,count); hold on; count = count + 1;
    fx = freq{i}.freq;
    py = squeeze(mean(freq{i}.powspctrm,1));
    pe = squeeze(sem(freq{i}.powspctrm));
    shadedErrorBar(fx,py(:,1),pe(:,1),{'color',col(1,:)},1);
    shadedErrorBar(fx,py(:,2),pe(:,2),{'color',col(2,:)},1);
    xlim([fx(1) fx(end)]); xlabel('Freq. (Hz)'); ylabel([lab{i},' Power (a.u.)'])
    title('Periodic Signal')

    subplot(3,4,count); hold on; count = count + 1;
    fx = freq{i}.freq;
    py = squeeze(mean(freq{i}.powspctrm(:,:,2) - freq{i}.powspctrm(:,:,1),1));
    pe = squeeze(sem(freq{i}.powspctrm(:,:,2) - freq{i}.powspctrm(:,:,1)));
    shadedErrorBar(fx,py,pe,{'color',col(3,:)},1);
    xlim([fx(1) fx(end)]); xlabel('Freq. (Hz)'); ylabel('Power Difference (a.u.)')
    title('Periodic Difference')

    subplot(3,4,count); hold on; count = count + 1;
    fy = squeeze(mean(freq{i}.fractal,1));
    fe = squeeze(sem(freq{i}.fractal));
    shadedErrorBar(fx,fy(:,1),fe(:,1),{'color',col(1,:)},1);
    shadedErrorBar(fx,fy(:,2),fe(:,2),{'color',col(2,:)},1);
    xlim([fx(1) fx(end)]); xlabel('Freq. (Hz)'); ylabel('Power (a.u.)')
    title('Aperiodic Signal')

    subplot(3,4,count); hold on; count = count + 1;
    fx = freq{i}.freq;
    fy = squeeze(mean(freq{i}.fractal(:,:,2) - freq{i}.fractal(:,:,1),1));
    fe = squeeze(sem(freq{i}.fractal(:,:,2) - freq{i}.fractal(:,:,1)));
    shadedErrorBar(fx,fy,fe,{'color',col(3,:)},1);
    xlim([fx(1) fx(end)]); xlabel('Freq. (Hz)'); ylabel('Power Difference (a.u.)')
    title('Aperiodic Difference')
end
end

function [stat,tbl] = int_runStats(data,parameters,tails)

% get key parameters
n_paras = numel(parameters);
n_datas = numel(data);

% predefine result cell
s = cell(n_paras,1);
t = cell(n_paras,1);

% cycle through each parameter
for i = 1 : n_paras

    % test parameter
    cfg             = [];
    cfg.tail        = zeros(n_datas,1)+tails(i);
    cfg.parameter   = parameters{i};
    [s{i},t{i}]     = run_oneSampleT(cfg, data(:));
end

% combine outputs
stat = cat(1,s{:});
tbl = cat(1,t{:});

% add labels
parameters = repmat(parameters,[n_datas 1]);
tbl = addvars(tbl,parameters(:),'Before','t');

end

function [data] = int_convertData(data,roi)

% get number of samples
n_subj = size(data.powspctrm,1);

% add in source model details
data.powdimord   = 'rpt_pos';
data.dim         = roi.dim;
data.inside      = roi.inside(:);
data.pos         = roi.pos;
data.pow         = zeros(n_subj,numel(data.inside));
data.slp         = zeros(n_subj,numel(data.inside));
data.cfg         = [];

% add in "outside" virtual electrodes to pow
tmp_pow = zeros(size(data.pow,1),size(data.inside,1));
tmp_pow(:,data.inside) = data.powspctrm;
data.pow = tmp_pow;

% add in "outside" virtual electrodes to slope
tmp_pow = zeros(size(data.pow,1),size(data.inside,1));
tmp_pow(:,data.inside) = data.slope;
data.slp = tmp_pow;

% remove freq details
data = rmfield(data,{'powspctrm','slope','label','freq','time','dimord','cfg'});
    
end

function tbl = int_extractSubjs(data,stat,tbl,metric,condition)

count = 1;
output = nan(size(data{1}.pow,1),numel(metric)*numel(condition));
label = cell(numel(metric)*numel(condition),1);
for i = 1 : numel(metric)
    for j = 1 : numel(condition)
        
        % define tail of cluster
        if tbl.t(count) > 0
            cluslab = 'posclusterslabelmat';
        elseif tbl.t(count) < 0
            cluslab = 'negclusterslabelmat'; 
        else
            label{count} = [condition{j},'_',metric{i}];
            count = count + 1;
            continue; 
        end
        
        % get cluster indices
        clusidx = stat{count}.(cluslab)==1;
        
        % extract data
        dat = nanmean(data{j}.(metric{i})(:,clusidx),2);
        
        % log normalise
        dat = dat ./ 10.^max(order(dat));
        sgn = sign(dat);
        logdat = log10(abs(dat)+1) .* sgn .* 10;
        output(:,count) = logdat;
        
        % define label
        label{count} = [condition{j},'_',metric{i}];
        count = count + 1;
    end
end

% create table
tbl = array2table(output,'VariableNames',label);

% plot
figure('visible','on'); hold on
subplot(1,2,1); hold on
boxplot(table2array(tbl(:,1:size(tbl,2)/2)))
for i = 1 : size(tbl,2)/2  
    plot(i+.3,table2array(tbl(:,i)),'ko');
end
subplot(1,2,2); hold on
boxplot(table2array(tbl(:,size(tbl,2)/2+1:end)))
for i = 1 : size(tbl,2)/2  
    plot(i+.3,table2array(tbl(:,i+size(tbl,2)/2)),'ko');
end

end

function int_extractMap(stat,filename)

% get indices of clusters
clus_idx = stat.negclusterslabelmat==1;

% create source data structure
source                  = [];
source.inside           = stat.inside;
source.dim              = stat.dim;
source.pos              = stat.pos*10;
source.unit             = 'mm';

% define powspctrm of cluster
source.pow              = nan(size(stat.pos,1),1);     
source.pow(clus_idx)	= stat.stat(clus_idx); 

% reshape data to 3D
source.pow              = reshape(source.pow,source.dim);
source.inside           = reshape(source.inside,source.dim);

% add transformation matrix
source.transform        = [1,0,0,-91;
                           0,1,0,-127;
                           0,0,1,-73;
                           0,0,0,1];

% export
cfg = [];
cfg.parameter     = 'pow';               % specify the functional parameter to write to the .nii file
cfg.filename      = filename;  % enter the desired file name
cfg.filetype      = 'nifti';
cfg.vmpversion    = 2;
cfg.datatype      = 'float';
ft_volumewrite(cfg,source);      % be sure to use your interpolated source data

% reslice to 1mm isometric to match template MRI
reslice_nii(filename,filename,[1 1 1]);

end
