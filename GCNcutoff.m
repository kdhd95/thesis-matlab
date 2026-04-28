function [gcn,A] = GCNcutoff(data,threshold)

corrtype = 'Spearman'; % the correlation type

C = corr(data','Type',corrtype); % the autocorrelation of the genes
C = abs(C - diag(diag(C))); % removing the diagonal - all ones.

A = C>=threshold;
links = sum(A,1);
gcn = mean(links);
end 