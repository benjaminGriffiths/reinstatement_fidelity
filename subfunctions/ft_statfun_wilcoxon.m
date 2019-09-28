function [s, cfg] = ft_statfun_wilcoxon(cfg, dat, design)

% set defaults
if ~isfield(cfg, 'computestat'),    cfg.computestat    = 'yes'; end
if ~isfield(cfg, 'computecritval'), cfg.computecritval = 'no';  end
if ~isfield(cfg, 'computeprob'),    cfg.computeprob    = 'no';  end
if ~isfield(cfg, 'alpha'),          cfg.alpha          = 0.05;  end
if ~isfield(cfg, 'tail'),           cfg.tail           = 1;     end

% perform some checks on the configuration
if strcmp(cfg.computeprob,'yes') && strcmp(cfg.computestat,'no')
    ft_error('P-values can only be calculated if the test statistics are calculated.');
end
if ~isfield(cfg,'uvar') || isempty(cfg.uvar)
    ft_error('uvar must be specified for dependent samples statistics');
end

% perform some checks on the design
sel1 = find(design(cfg.ivar,:)==1);
sel2 = find(design(cfg.ivar,:)==2);
n1  = length(sel1);
n2  = length(sel2);
if (n1+n2)<size(design,2) || (n1~=n2)
    ft_error('Invalid specification of the design array.');
end
nunits = length(design(cfg.uvar, sel1));
df = nunits - 1;
if nunits<2
    ft_error('The data must contain at least two units (usually subjects).')
end
if (nunits*2)~=(n1+n2)
    ft_error('Invalid specification of the design array.');
end

if strcmp(cfg.computestat,'yes')
    % compute the statistic
    % store the positions of the 1-labels and the 2-labels in a nunits-by-2 array
    poslabelsperunit = zeros(nunits,2);
    poslabel1        = find(design(cfg.ivar,:)==1);
    poslabel2        = find(design(cfg.ivar,:)==2);
    [~,i]          = sort(design(cfg.uvar,poslabel1), 'ascend');
    poslabelsperunit(:,1) = poslabel1(i);
    [~,i]          = sort(design(cfg.uvar,poslabel2), 'ascend');
    poslabelsperunit(:,2) = poslabel2(i);
    
    % calculate the differences between the conditions
    diffmat = dat(:,poslabelsperunit(:,1)) - dat(:,poslabelsperunit(:,2));
    
    % cycle through channels
    for chan = 1 : size(diffmat,1)
    
        % get sign
        diffvec = diffmat(chan,:);
        signvec = diffvec;
        signvec(signvec>=0) = 1;
        signvec(signvec<0)  = -1;
        
        % sort based on absolute difference
        [~,idx] = sort(abs(diffvec));
        
        % calculate wilcoxon sign
        W = sum(signvec(idx) .* idx);
        
        % convert to zscore
        n = numel(signvec);
        oW = sqrt((n*(n+1)*(2*n+1))/6);
        s.stat(chan,1) = W ./ oW;
    end
end

if strcmp(cfg.computecritval, 'yes')
    % also compute the critical values
    s.df      = df;
    if cfg.tail==-1
        s.critval = tinv(cfg.alpha,df);
    elseif  cfg.tail==0
        s.critval = [tinv(cfg.alpha/2,df),tinv(1-cfg.alpha/2,df)];
    elseif cfg.tail==1
        s.critval = tinv(1-cfg.alpha,df);
    end
end

if strcmp(cfg.computeprob, 'yes')
    % also compute the p-values
    s.df      = df;
    if cfg.tail==-1
        s.prob = tcdf(s.stat,s.df);
    elseif  cfg.tail==0
        s.prob = 2*tcdf(-abs(s.stat),s.df);
    elseif cfg.tail==1
        s.prob = 1-tcdf(s.stat,s.df);
    end
end