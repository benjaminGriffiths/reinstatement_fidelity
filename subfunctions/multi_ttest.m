function [pval,stat] = multi_ttest(X)

[n,p]   = size(X);
mu      = zeros(1,p);
m       = mean(X); % Mean vector from data matrix X.
S       = cov(X);  % Covariance matrix from data matrix X.
T2      = n*(m-mu)*inv(S)*(m-mu)'; %Hotelling's T-Squared statistic.

% run F approximation
F       = (n-p)/((n-1)*p)*T2;
v1      = p;  %Numerator degrees of freedom.
v2      = n-p;  %Denominator degrees of freedom.
pval    = 1 - fcdf(F,v1,v2);  %Probability that null Ho: is true.

% create output structure
stat.F = F;
stat.v1 = v1;
stat.v2 = v2;
stat.T2 = T2;