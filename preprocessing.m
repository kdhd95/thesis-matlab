function [type,young,old] = preprocessing(filey,path,numofgenes)
% This function is performing the import and preprocessing to the datasets

% Importing - the user can only pick the young datasets
% [filey,path] = uigetfile({'*Young*.mat'},'choose young data','/Users/Daniel1/תואר שני/תזה/data/');

young = load([path filey]);
young = young.data2;

% Changing the file name for finding the coresponding old dataset
fileo = replace(filey,"young","old");
fileo = replace(fileo,"Young","Old");
% Importing the old dataset
old = load([path fileo]);
old = old.data1;
type = filey(12:end-4);


% Checking if the number of cells is matching bitween the datasets
numcelly = size(young,2);
numcello = size(old,2);
if numcelly~=numcello
    %finding the low number of cells in datasets
    minum = min(numcelly,numcello);
    % throw out the number of extra cells in random
    old = old(:,randperm(numcello,minum));
    young = young(:,randperm(numcelly,minum));
end

% Reducing the number of genes as requested
if numofgenes~=0
    [young,old] = reduce(young,old,numofgenes);
end


