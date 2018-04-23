function [ freeSpacePathLosses, losLengths ] ...
    = computeFreeSpacePathLosses(tx3D, rxs3D, fCarrierGHz)
%COMPUTEFREESPACEPATHLOSSES Compute the free-space path losses in dB
%between one TX location and many RX locations.
%
% Input:
%   - tx3D
%     A three-element vector specifying the locations [x y altitude] of the
%     transmitter
%   - rxs3D
%     A matrix with each row representing a RX location in the form of [x y
%     altitude].
%   - fCarrierGHz
%     The signal frequency in GHz.
% Outputs:
%   - freeSpacePathLosses
%     A column vector holding the free-space path losses in dB, which each
%     element cooresponding to one RX location in rxs3D.
%   - losLengths
%     A column vector holding the TX-to-RX distances, which each element
%     cooresponding to one RX location in rxs3D.
%
% Yaguang Zhang, Purdue, 04/20/2018

lambda = physconst( 'LightSpeed' )/(fCarrierGHz*10^9);

[numRxLocs, ~] = size(rxs3D);
losLengths = vecnorm((rxs3D - repmat(tx3D, numRxLocs, 1))')';

freeSpacePathLosses = fspl(losLengths, lambda);

end

% EOF