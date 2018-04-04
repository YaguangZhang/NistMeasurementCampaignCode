function [ hpbw, anglesHpbw, powerInDbHpbw ] ...
    = computeHPBW( angles, powersInDb, maxPowerInDb )
%COMPUTEHPBW Compute the Half Power Beamwidth for an antenna patter.
%
% Inputs:
%   - angles, powersInDb
%     Column vectors specifying the pattern of the antenna. The angles
%     vector should start with 0 and end with 360 in degree. The powersInDb
%     vector contains powers in dB.
%   - maxPowerInDb
%     Optional. A scalar for the expected maximum power of the pattern in
%     dB.
%
% Outputs:
%   - hpbw
%     The resulted Half Power Beamwidth in degree.
%   - anglesHpbw, powerInDbHpbw
%     The two power pattern points which coorespond to the HPBW.
%
% Yaguang Zhang, Purdue, 10/12/2017

if nargin < 3
    maxPowerInDb = max(powersInDb);
end
anglesHpbw = nan(2,1);

% Harf Power => -3 dB.
powerInDbShreshold = maxPowerInDb-3;

% Extend the angles range with this factor. We will linearly fit more power
% points for computing the results.
RATIO_NUM_OF_ANGLES = 5;

assert(angles(1)==0 && angles(end)==360, ...
    'The angles vector should start with 0 and end with 360 in degree!');
assert(powersInDb(1)>powerInDbShreshold && angles(end)>powerInDbShreshold, ...
    'The powers vector should both start with and end with a value larger than powerInDbShreshold!');

extAngles = (0:(360/(length(angles)-1)/RATIO_NUM_OF_ANGLES):360)';
extPowersInDb = interp1q(angles, powersInDb, extAngles);

indicesHpbw = [find(extPowersInDb<=powerInDbShreshold, 1, 'first'); ...
    find(extPowersInDb<=powerInDbShreshold, 1, 'last')];
anglesHpbw = extAngles(indicesHpbw);
powerInDbHpbw = extPowersInDb(indicesHpbw);

hpbw = anglesHpbw(1) + 360 - anglesHpbw(2);
end
%EOF