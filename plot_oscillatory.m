%% Oscillatory Spectrum
% figure
figure; 
hold on

% plot parameters
freq = group_irasa{subj,i,1}.freq;
col = [0.8 0.4 0.4; 0.4 0.4 0.8];

% cycle through subjects
for subj = 1 : n_subj

    % subplot
    subplot(6,4,subj);
    hold on
    
    % cycle through encoding and retrieval
    for i = 1 : 2
    
        % get average oscillatory spectrum
        osci_spec = cat(1,group_irasa{subj,i,1}.prepow,group_irasa{subj,i,1}.powspctrm);
        osci_spec = squeeze(mean(osci_spec,1));
        
        % plot
        shadedErrorBar(freq,mean(osci_spec),sem(osci_spec),{'color',col(i,:)},1);
    end
    
%     % get average oscillatory spectrum
%     osci_spec = cat(1,group_irasa{subj,1,1}.prepow,group_irasa{subj,1,1}.powspctrm,group_irasa{subj,2,1}.prepow,group_irasa{subj,2,1}.powspctrm);
%     osci_spec = squeeze(mean(osci_spec,1));
% 
%     % plot
%     shadedErrorBar(freq,mean(osci_spec),sem(osci_spec),[],1);

    % figure details
    xlim([freq(1) freq(end)])
    set(gca,'xtick',linspace(freq(1),freq(end),5),'xticklabel',5:5:25)
    title(sprintf('sub-%02.0f [red is perception, blue is retrieval]',subj));
end

%% Peak Histogram
% figure
ax = tight_subplot(7,3,[0.05 0.04]);
hold on

% plot parameters
freq = group_irasa{subj,i,1}.freq;
col = [0.8 0.4 0.4; 0.4 0.4 0.8];
            
% cycle through subjects
for subj = 1 : n_subj

    % subplot
    peak_freq = cell(2);
    axes(ax(subj))
    
    % cycle through encoding and retrieval
    for i = 1 : 2

        % cycle through channels
        for chan = 1 : n_chan

            % get average oscillatory spectrum
            osci_spec = cat(1,group_irasa{subj,i,1}.prepow(:,chan,:),group_irasa{subj,i,1}.powspctrm(:,chan,:));
            osci_spec = squeeze(mean(osci_spec,1));
            
            % find peaks of spectrum
            [vals,idx] = findpeaks(osci_spec);
            
            % if peaks are present
            if ~isempty(idx)
                [~,midx] = max(vals);
                idx = idx(midx);

                % add peak power to regressor table    
                peak_freq{i}(chan) = freq(idx);
            end            
        end     
        
        % plot      
        histogram(peak_freq{i},5:25,'normalization','probability');
        hold on
    end
    
    % plot details    
    set(ax(subj),'box','off','tickdir','out','xtick',6.5:4:24.5,'xticklabel',6:4:24)
    title(ax(subj),sprintf('sub-%02.0f',subj));
    if ismember(subj,[19 20 21]); xlabel(ax(subj),'Freq. (Hz)'); end
    if ismember(subj,[1 4 7 10 13 16 19]); ylabel(ax(subj),'Probability'); end
    ylim(ax(subj),[0 1])
end


%% Peak Histogram
% figure
ax = tight_subplot(7,3,[0.05 0.04]);
hold on

% plot parameters
freq = group_irasa{subj,i,1}.freq;
col = [0.8 0.4 0.4; 0.4 0.4 0.8];
            
% cycle through subjects
for subj = 1 : n_subj

    % subplot
    peak_freq = cell(2);
    axes(ax(subj))
    
    % cycle through channels
    for chan = 1 : n_chan

        % get average oscillatory spectrum
        osci_spec = cat(1,group_irasa{subj,1,1}.prepow(:,chan,:),group_irasa{subj,1,1}.powspctrm(:,chan,:),group_irasa{subj,2,1}.prepow(:,chan,:),group_irasa{subj,2,1}.powspctrm(:,chan,:));
        osci_spec = squeeze(mean(osci_spec,1));

        % find peaks of spectrum
        [vals,idx] = findpeaks(osci_spec);

        % if peaks are present
        if ~isempty(idx)
            [~,midx] = max(vals);
            idx = idx(midx);

            % add peak power to regressor table    
            peak_freq{i}(chan) = freq(idx);
        end            
    end     

    % plot      
    histogram(peak_freq{i},5:25,'normalization','probability','facecolor',[0.5 0.5 0.5]);
    hold on
    
    % plot details    
    set(ax(subj),'box','off','tickdir','out','xtick',6.5:4:24.5,'xticklabel',6:4:24)
    title(ax(subj),sprintf('sub-%02.0f',subj));
    if ismember(subj,[19 20 21]); xlabel(ax(subj),'Freq. (Hz)'); end
    if ismember(subj,[1 4 7 10 13 16 19]); ylabel(ax(subj),'Probability'); end
    ylim(ax(subj),[0 1])
end
