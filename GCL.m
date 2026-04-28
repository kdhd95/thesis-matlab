function [gcl] = GCL(X)
% the GCL function that get the data <X>, and calculating the GCL with
% bcdCorr method as the correlation.
% 

N = size(X,1);
m = 50;
D = zeros(m,1);

for i=1:m
    idx = randperm(N);
    S1 = X(idx(1:floor(N/2)),:)';
    S2 = X(idx((floor(N/2)+1):end),:)';
    D(i) = dCov(S1,S2)./sqrt(dCov(S1,S1).*dCov(S2,S2));
end

gcl = sum(D)./m;
end

