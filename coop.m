function coop_sum = coop(x, w_coop)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = length(x);
coop_sum = zeros(N, 1);
parfor j = 1:N
    for k = (j+1):N  % Ensure j < k
        coop_sum = coop_sum + ...
            w_coop(:,j,k).*(x(j)./(1+x(j))).*(x(k)./(1+x(k)));
    end
end    
end
