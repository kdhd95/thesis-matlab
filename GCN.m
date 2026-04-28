function [gcn,A,C] = GCN(data,alpha)
% This function is calaulating the GCN method as I implement it. 
% the input is the data.
% The output is the gcn index.

% alpha = 10^-4;
corrtype = 'Spearman'; % the correlation type
Mu = zeros(size(data,1));
M2 = zeros(size(data,1));

C = corr(data','Type',corrtype); % the autocorrelation of the genes
C = abs(C - diag(diag(C))); % removing the diagonal - all ones.

for i=1:100
    random = randomized(data);
    Cr = corr(random','Type',corrtype); % the autocorrelation of the genes
    % Welford's online algorithm
    delta = Cr - Mu;
    Mu = Mu + delta./i;
    delta2 = Cr - Mu;
    M2 = M2 + delta .* delta2;
end

Sigma =sqrt(M2./99);
P_values = 1-normcdf(C,Mu,Sigma);
A = P_values<=alpha;
links = sum(A,1);
gcn = mean(links);
end