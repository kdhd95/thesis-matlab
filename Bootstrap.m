function [bootstat] = Bootstrap(funcr,M,numsub,persub)
% This function doing bootstrap for a function <funcr> with data <M>,
% doing it <numsub> times, each time with <persub> percent from M.

numcells = size(M,2); % number of cells for to choose from.
bootstat = zeros(1,numsub); % the vector of results

parfor i=1:numsub
    index = randperm(numcells); % randomly picking the cells.
    data = M(:,index(1:floor(numcells*persub))); % lessen the data.
    bootstat(i) = funcr(data); % perform the function
end

end