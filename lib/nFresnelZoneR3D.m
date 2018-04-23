function [ nRadius ] = nFresnelZoneR3D(n, lambda, tx3D, rx3D, p3D)
%NFRESNELZONER3D Compute the radius of the n-th Fresenl zone at a given
%point in 3D.
%
% Input:
%   - n
%     The n-th Fresnel zone to consider.
%   - lambda
%     The wavelength of the transmitted signal.
%   - tx3D, rx3D, p3D
%     Three-element vectors specifying the locations [x y z] of the
%     transmitter, reciever and the observation point, respectively.
% Output:
%   - nRadius
%     The radius of the n-th Fresenl zone at the given point in 3D. Note
%     that it will be 0 for any point that is not between the endpoints of
%     the link.
%
% Note: all the length parameters should have the same unit.
%
% Ref: https://en.wikipedia.org/wiki/Fresnel_zone
%
% Yaguang Zhang, Purdue, 04/18/2018

% Make sure the input 3D points are row vectors.
assert(isrow([tx3D rx3D p3D]), ...
    'The input 3D points should all be row vectors!');

% LOS line.
los = createLine3d(tx3D, rx3D);

% The position of the observation point onto the LOS line.
pOnLos3D = los(1:3) + linePosition3d(p3D, los) .* los(4:6);

% Only output positive radius if p is between the transmitter and the
% receiver.
deltaTxToPOnLos = pOnLos3D-tx3D;
deltaPToRxOnLos = rx3D-pOnLos3D;
if all(sign(deltaTxToPOnLos) == sign(deltaPToRxOnLos))
    % Distance between TX and pOnLos3D.
    d1 = norm(deltaTxToPOnLos);
    % Distance between pOnLos3D and RX.
    d2 = norm(deltaPToRxOnLos);
    nRadius = sqrt( (n * lambda * d1 * d2) / (d1 + d2) );
else
    nRadius = 0;
end

end
% EOF