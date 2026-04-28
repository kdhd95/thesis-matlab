function [RandData] = randomized(Data)
% This function is randomizing a data set in the genes way.
% after this random, there should be no connections between the genes in a
% single cell.

[genesNumber,cellNumber] = size(Data); 
RandData = zeros(genesNumber,cellNumber); %the new random data

for r = 1:genesNumber
    % randomizing each row so the genes will mix bitween different cells
    RandData(r,:) = Data(r,randperm(cellNumber)); 
end

end
