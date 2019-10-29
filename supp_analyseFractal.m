function supp_analyseFractal

% define samples
vals = [0.1:0.1:0.9,1:10,20:10:100];

% predefine fit
fracFit = zeros(numel(vals),100,2);

% cycle through each sample
for samp = 1 : numel(vals)

    % load data
    dat = csvread(['E:\bjg335\projects\reinstatement_fidelity\data\supp_fractal\sig_',num2str(samp),'.csv']);
    
    % get fractal
    output = amri_sig_fractal(dat,100);
    
    % cycle through reps
    for trl = 1 : size(output.frac,2)

        % log transform fractal
        logF = log10(output.freq);
        logP = log10(output.frac(:,trl));

        % get linfit
        tmp = fitlm(logF,logP); %post

        % get slope and intercept difference
        fracFit(samp,trl,:) = tmp.Coefficients.tStat;
    end
    
    % updare
    fprintf('part %02.0f of %02.0f complete...\n',samp,numel(vals))
end

% calculate signal to noise
snr = squeeze(abs(mean(fracFit,2) ./ std(fracFit,[],2)));

% plot
figure; hold on
plot(vals,snr);
set(gca,'xscale','log')
xlabel('Epoch Duration (s)')
ylabel('Signal-to-Noise Ratio')
legend('intercept','slope','location','northwest')