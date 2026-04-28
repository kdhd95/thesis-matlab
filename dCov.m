function [d] = dCov(X,Y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

M = size(X,1);

a = pdist2(X,X);
al = sum(a,2)./M;
ar = sum(a,1)./M;
aa = sum(a,'all')./M^2;

b = pdist2(Y,Y);
bl = sum(b,2)./M;
br = sum(b,1)./M;
bb = sum(b,'all')./M^2;

A = (a-al-ar+aa-a./M).*M/(M-1);
B = (b-bl-br+bb-b./M).*M/(M-1);

for i=1:M
    A(i,i) = (al(i)-aa).*M/(M-1);
    B(i,i) = (bl(i)-bb).*M/(M-1);
end

d = (sum(A.*B,'all')-trace(A.*B).*M/(M-2))./(M*(M-3));


end

