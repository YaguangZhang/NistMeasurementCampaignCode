function [distFromTreeToTx, distFromTreeToRx] ...
    = distsFromTreeToTerminals3D(tree3D, tx3D, rx3D)
%DISTSFROMTREETOTERMINALS3D Estimates the distances from a tree to the TX
%and the RX in 3D.
%   Inputs:
%       - tree3D, tx3D, rx3D
%         The 3D locations for the tree, tx and the rx in the form of [x y
%         altitude].
%   Outputs:
%       - distFromTreeToTx, distFromTreeToRx
%         The distances from the tree to the TX and the RX, respectively.
%
% Yaguang Zhang, Purdue, 09/17/2018

computeDist = @(x,y) sqrt(sum((x-y).^2));

distFromTxToRx3D = computeDist(rx3D,tx3D);

distFromTreeToTx2D = computeDist(tree3D(1:2), tx3D(1:2));
distFromTreeToRx2D = computeDist(tree3D(1:2), rx3D(1:2));

ratio2dDistTreeToTx = distFromTreeToTx2D ...
    ./(distFromTreeToTx2D+distFromTreeToRx2D);

distFromTreeToTx = distFromTxToRx3D.*ratio2dDistTreeToTx;
distFromTreeToRx = distFromTxToRx3D-distFromTreeToTx;

end

% EOF