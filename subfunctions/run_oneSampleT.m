function [stat,tbl] = run_oneSampleT(cfg,data)

% set default parameters if not input
if ~isfield(cfg,'frequency'); cfg.frequency = [8 30]; end
if ~isfield(cfg,'latency'); cfg.latency = [0.5 1.5]; end
if ~isfield(cfg,'parameter'); cfg.parameter = 'powspctrm'; end
if ~isfield(cfg,'statistic'); cfg.statistic = 'ft_statfun_depsamplesT'; end
if ~isfield(cfg,'tail'); cfg.tail = 0; cfg.clustertail = 0; end
if ~isfield(cfg,'rm_outliers'); cfg.rm_outliers = false; end

% if data is struct, convert to cell
if isstruct(data); data = {data}; end

% get number of data structures to analyse
n_data = numel(data);

% predefine cell for statistics
stat = cell(n_data,1);

% prepare table for stat values
tbl = array2table(nan(n_data,4),'VariableNames',{'t','p','ci','d'});

% define statisitics
config                   = [];
config.uvar              = 1;
config.ivar              = 2;
config.method            = 'montecarlo';
config.statistic         = cfg.statistic;
config.correctm          = 'cluster';
config.numrandomization  = 2000;
config.alpha             = 0.05;
config.correcttail       = 'prob';
config.parameter         = 'pow';

% cycle through data structures
for condition = 1 : n_data

    % create pow parameter
    data{condition}.pow = data{condition}.(cfg.parameter);
    
    % define tail    
    config.tail = cfg.tail(condition);
    
    % predefine conditionals
    issource = false;
    issingle = false;
            
    % check if there is only a single comparison
    if size(data{condition}.(cfg.parameter),2) == 1 && size(data{condition}.(cfg.parameter),3) == 1 && size(data{condition}.(cfg.parameter),4) == 1
        issingle = true;
        fprintf('Assuming data has a single comparison...\n')
       
    % check whether data is in source space
    elseif ~isfield(data{condition},'label')
        fprintf('No label information detected, assuming data is in source space...\n')
        issource = true;
        
    else
        fprintf('Assuming data is channel level...\n')
    end
    
    % define statistical design
    n_subj      = size(data{condition}.(config.parameter),1);
    config.design  = [1:n_subj, 1:n_subj; ones(1,n_subj), ones(1,n_subj)+1];
       
    % create null hypothesis
    null_freq = create_null_hypothesis(data{condition},config.parameter);
    
    % if data is in source space
    if issource && ~issingle
        
        % run statistics
        config.dim          = data{condition}.dim;  % specify dimensions of your source grid        
        config.clusteralpha = 0.05;
        config.clustertail  = cfg.tail(condition);
        stat{condition}     = ft_sourcestatistics(config, data{condition}, null_freq);
        
    % if data consists of a single comparison
    elseif issingle
        
        % run statistics
        config.correctm     = 'no'; % remove multiple comparisons
        stat{condition}     = ft_freqstatistics(config, data{condition}, null_freq);
        
    else
        
        % prepare neighbours  
        tmp                 = [];
        tmp.layout          = cfg.layout;
        tmp.method          = 'triangulation';
        config.neighbours   = ft_prepare_neighbours(tmp);

        % run statistics
        config.frequency    = cfg.frequency;
        config.latency      = cfg.latency;        
        config.clusteralpha = 0.05;
        config.clustertail  = cfg.tail;
        stat{condition}     = ft_freqstatistics(config, data{condition}, null_freq);
    end
     
    % for coding help, add in missing fields
    if ~isfield(stat{condition},'posclusters'); stat{condition}.posclusters = []; end
    if ~isfield(stat{condition},'negclusters'); stat{condition}.negclusters = []; end
    
    % if one-tailed test
    if cfg.tail(condition)~=0

        % define clusters of interest
        if cfg.tail(condition)==1
            tailname = 'posclusters';
        elseif cfg.tail(condition)==-1
            tailname = 'negclusters';
        end

        % if data does not consist of a single value
        if ~issingle

            % check cluster exists
            if isfield(stat{condition,1},tailname) && ~isempty(stat{condition,1}.(tailname))
                
                % extract key values
                tbl.p(condition,1)  = round(stat{condition,1}.(tailname)(1).prob,3);
                tbl.t(condition,1)  = round(stat{condition,1}.(tailname)(1).clusterstat ./ sum(stat{condition,1}.([(tailname),'labelmat'])(:) == 1),3);
                tbl.ci(condition,1) = round(stat{condition,1}.(tailname)(1).cirange,3);

                % calculate cohen's dz
                tbl.d(condition,1) = round(tbl.t(condition,1)./ sqrt(n_subj),3);
            end

        else
            % extract key values
            tbl.p(condition,1)  = round(stat{condition,1}.prob,3);
            tbl.t(condition,1)  = round(stat{condition,1}.stat,3);
            tbl.ci(condition,1) = round(stat{condition,1}.cirange,3);

            % if parametric, calculate cohen's dz
            tbl.d(condition,1) = round(tbl.t(condition,1)./ sqrt(n_subj),3);
        end
        
    else
        % if data does not consist of a single value
        if ~issingle % CURRENTLY NOT FUNCTIONAL

            
            
            % identify which p-value is larger
            if isempty(stat{condition,1}.posclusters) && ~isempty(stat{condition,1}.negclusters)
                tailname = 'negclusters';
            elseif ~isempty(stat{condition,1}.posclusters) && isempty(stat{condition,1}.negclusters)
                tailname = 'posclusters';
            elseif ~isempty(stat{condition,1}.posclusters) && ~isempty(stat{condition,1}.negclusters)
                if stat{condition,1}.posclusters(1).prob < stat{condition,1}.negclusters(1).prob
                    tailname = 'posclusters';
                else                    
                    tailname = 'negclusters';
                end
            else
                continue
            end
            
            % extract key values
            tbl.p(condition,1)  = round(stat{condition,1}.(tailname)(1).prob,3);
            tbl.t(condition,1)  = round(stat{condition,1}.(tailname)(1).clusterstat ./ sum(stat{condition,1}.([(tailname),'labelmat'])(:) == 1),3);
            tbl.ci(condition,1) = round(stat{condition,1}.(tailname)(1).cirange,3);

            % calculate cohen's dz
            tbl.d(condition,1) = round(tbl.t(condition,1)./ sqrt(n_subj),3);

        else
            % extract key values
            tbl.p(condition,1)  = round(stat{condition,1}.prob,3);
            tbl.t(condition,1)  = round(stat{condition,1}.stat,3);
            tbl.ci(condition,1) = round(stat{condition,1}.cirange,3);

            % calculate cohen's dz
            tbl.d(condition,1) = round(tbl.t(condition,1)./ sqrt(n_subj),3);
        end
    end
end
