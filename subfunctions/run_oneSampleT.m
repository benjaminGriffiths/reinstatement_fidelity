function [stat,tbl] = run_oneSampleT(cfg,data)

% set default parameters if not input
if ~isfield(cfg,'frequency'); cfg.frequency = [8 30]; end
if ~isfield(cfg,'latency'); cfg.latency = [0.5 1.5]; end
if ~isfield(cfg,'parameter'); cfg.parameter = 'powspctrm'; end
if ~isfield(cfg,'tail'); cfg.tail = 0; cfg.clustertail = 0; end

% if data is struct, convert to cell
if isstruct(data); data = {data}; end

% get number of data structures to analyse
n_data = numel(data);

% predefine cell for statistics
stat = cell(n_data,1);

% prepare table for stat values
tbl = array2table(zeros(n_data,4),'VariableNames',{'t','p','ci','d'});

% define statisitics
config                   = [];
config.uvar              = 1;
config.ivar              = 2;
config.method            = 'montecarlo';
config.statistic         = 'ft_statfun_depsamplesT';
config.correctm          = 'cluster';
config.clusteralpha      = 0.05;
config.numrandomization  = 2000;
config.alpha             = 0.05;
config.tail              = cfg.tail;
config.clustertail       = cfg.tail;
config.parameter         = cfg.parameter;

% cycle through data structures
for condition = 1 : n_data

    % define statistics type
    if ~isfield(data{condition},'label')
        fprintf('No label information detected, assuming data is in source space...\n')
        issource = true;
    else
        fprintf('Assuming data is channel level...\n')
        issource = false;
    end
    
    % define statistical design
    n_subj      = size(data{condition}.(config.parameter),1);
    config.design  = [1:n_subj, 1:n_subj; ones(1,n_subj), ones(1,n_subj)+1];
       
    % create null hypothesis
    null_freq = create_null_hypothesis(data{condition},config.parameter);
    
    % if data is in source space
    if issource
        
        % run statistics
        config.dim         = data{condition}.dim;  % specify dimensions of your source grid
        stat{condition} = ft_sourcestatistics(config, data{condition}, null_freq);
    
    else
        
        % prepare neighbours  
        config             = [];
        config.layout      = lay;
        config.method      = 'triangulation';
        config.neighbours  = ft_prepare_neighbours(config);

        % run statistics
        config.frequency   = cfg.frequency;
        config.latency     = cfg.latency;
        stat{condition}     = ft_freqstatistics(config, data{condition}, null_freq);
    end
        
    % extract key values
    tbl.p(condition,1)  = round(stat{condition,1}.negclusters(1).prob,3);
    tbl.t(condition,1)  = round(stat{condition,1}.negclusters(1).clusterstat ./ sum(stat{condition,1}.negclusterslabelmat(:) == 1),3);
    tbl.ci(condition,1) = round(stat{condition,1}.negclusters(1).cirange,3);

    % calculate cohen's dz
    tbl.d(condition,1) = round(tbl.t(condition,1)./ sqrt(n_subj),3);
end
