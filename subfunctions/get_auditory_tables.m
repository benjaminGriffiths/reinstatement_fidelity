function [stim_details,scan_details] = get_auditory_tables(dir_root,subj_handle)

% define key parameters
n_trials = 192;
n_volumes = 255;
n_runs = 8;

% create table to record stimulus detail
stim_details  = array2table(zeros(n_trials*2,3),'VariableNames', {'encoding','modality','memory'});
stim_count = 1;

% create table to record scan details
scan_details  = array2table(zeros(n_volumes*8,2),'VariableNames', {'encoding','modality'});
scan_count = 1;  

% cycle through each run
for run = 1 : n_runs

    % load event table
    tbl = readtable([dir_root,subj_handle,'/func/',subj_handle,'_task-rf_run-',num2str(run),'_events.tsv'],'FileType','text','Delimiter','\t');

    % check block types
    block_encoding  = double(any(ismember(tbl.operation,'encoding')));
    block_visual    = double(any(ismember(tbl.modality,'Auditory')));

    % cycle through every event
    for e = 1 : size(tbl,1)

        % if an event
        if strcmpi(tbl.trial_type(e),'Stimulus Onset')

            % add key values
            stim_details.encoding(stim_count) = strcmpi(tbl.operation(e),'encoding');
            stim_details.modality(stim_count) = strcmpi(tbl.modality(e),'Auditory');

            % get stimulus values
            switch tbl.stimulus{e}
                case 'TRUMPET';     stim_details.stimulus(stim_count,1) = 1;
                case 'PIANO';       stim_details.stimulus(stim_count,1) = 2;
                case 'GUITAR';      stim_details.stimulus(stim_count,1) = 3;
                case 'ACCORDIAN';   stim_details.stimulus(stim_count,1) = 4;
            end

            % change value if retrieval
            if stim_details.encoding(stim_count) ~= 1 && stim_details.modality(stim_count) == 1
                stim_details.stimulus(stim_count,1) = stim_details.stimulus(stim_count,1) + 4;
            end

            % get memory
            stim_details.memory(stim_count) = tbl.recalled(e);
            
            % update counter
            stim_count = stim_count + 1; 

        % else if a volume
        elseif strncmpi(tbl.trial_type(e),'Volume',6)

            % add key values
            scan_details.encoding(scan_count) = block_encoding;
            scan_details.modality(scan_count) = block_visual;

            % update counter
            scan_count = scan_count + 1;                 
        end          
    end
end

% kick out first three/last five scans of each run
scan_details([1:3 251:258 506:513 761:768 1016:1023 1271:1278 1526:1533 1781:1788 2036:2040],:) = [];

% clean up
clear run tbl e stim_count
