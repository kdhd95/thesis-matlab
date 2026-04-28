function [young,old] = reduce(young,old,num)
% This function is reducing 2 datasets without ruin their compatibility
% between the genes. We are keeping the top <num> genes that are most
% present, with 5% error rate.

if (num*1.05)<size(young,1) %for checking if even necesarry
    
    % getting the present values of each gene.
    avy = sum(young,2);
    avo = sum(old,2);
    
    % sorting the genes for the most present ones.
    [~,Iy] = sort(avy,'descend');
    [~,Io] = sort(avo,'descend');

    % for doing the reduce for both datasets, we will compare the genes
    % that are most present between the datasets, and we will only take the
    % genes that are most present in both datasets.
    % For exaple: if i want only 100 genes. and gene A is ranked 50 in the 
    % young dataset, and 150 in the old dataset, than he wont make the
    % final cut. insted, gene 101 in the young dataset that is ranked 99 in
    % the old dataset will make the cut.

    numtry = num;
    cut = 0;
    while cut<(num) % we need to achive at least <num> genes
        %each loop we checking more genes to fulfill the requirment of <num>
        numtry = numtry+num*0.05;
        % taking the top genes
        tryy = Iy(1:numtry);
        tryo = Io(1:numtry);
        % checking that the genes are the same in both datasets
        indexs = intersect(tryy,tryo);
        cut = length(indexs); % number of aprroved genes
    end
    
    % after we decide what to throw, the reducing itself
    young = young(indexs,:);
    old = old(indexs,:);

end
end

