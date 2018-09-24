function [weightedTrunkNum, ...
    curDistsFromTreeToTx, curDistsFromTreeToRx, curNumTree3D] ...
    = computeWeightedTrunkNum( ...
    trees3D, txs3D, rx3D, weightFct)
%COMPUTEWEIGHTEDTRUNKNUM Compute the weighted trunk number.
%   We will compute the number of trunks in a weighted manner.
%
%   Inputs:
%       - trees3D, txs3D
%         Matrices representing the tree locations and the recevier
%         locations, respectively. Each row is one location in the form of
%         [x, y, altitude].
%       - rx3D
%         A row vector for the location of the receiver in the form of [x,
%         y, altitude].
%
%   Output:
%       - weightedTrunkNums
%         The resulting weighted trunk number:
%           sumForAllTrees(1/distToTx+1/distToRx)
%
% Yaguang Zhang, Purdue, 09/18/2018

[numRxLocs, ~] = size(txs3D);
weightedTrunkNum = 0;

if ~exist('weightFct', 'var')
    weightFct = @(d) -log(atan(0.1./d));
end

[curNumTree3D, ~]= size(trees3D);
[curDistsFromTreeToTx, curDistsFromTreeToRx] ...
    = deal(nan(curNumTree3D, 1));

if ~isempty(trees3D)    
    for idxRx = 1:numRxLocs
        [curDistFromTreeToTx, curDistFromTreeToRx] ...
            = distsFromTreeToTerminals3D(trees3D, txs3D(idxRx, :), rx3D);
        weightedTrunkNum = weightedTrunkNum ...
            + weightFct(curDistFromTreeToTx) ...
            + weightFct(curDistFromTreeToRx);
        
        curDistsFromTreeToTx(idxRx) = curDistFromTreeToTx;
        curDistsFromTreeToRx(idxRx) = curDistFromTreeToRx;
    end
end

end

% EOF