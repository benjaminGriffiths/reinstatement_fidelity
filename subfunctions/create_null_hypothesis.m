function null_freq = create_null_hypothesis(freq,parameter)

if nargin == 1
    parameter = 'powspctrm';
end

null_freq = freq;
null_freq.(parameter) = zeros(size(null_freq.(parameter)));