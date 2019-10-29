function supp_correlationVariability

%% Run Simulation
noise = 0 : 0.01 : 1;
ssivar = 0 : 0.01 : 1;
ntrls = 100;
r = nan(numel(noise),numel(ssivar));

for i = 1 : numel(noise)
    for j = 1 : numel(ssivar)
        
        % set ssi as ones plus variance
        ssi = ones(ntrls,1) + (rand(ntrls,1)*ssivar(j));
        
        % set eeg as ssi plus noise
        eeg = ssi + (rand(ntrls,1)*noise(i));
        
        % correlate
        r(i,j) = corr(ssi,eeg);
    end
end

% plot
imagesc(ssivar(2:end),noise,r(:,2:end)); axis xy
xlabel('SSI Variability'); ylabel('EEG Noise')

