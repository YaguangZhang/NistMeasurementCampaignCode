function [ beamLonLatPolyshape] = genBeamPolyshape(txLon, txLat, ...
    txAzInDeg, beamwidthInDeg, beamShapeSideLengthInM)
%GENBEAMPOLYSHAPE Generate a triangle polyshape to represent the area
%covered by the TX antenna beam.
%
% Inputs: 
%   - txLat, txLon, txAzInDeg
%     The GPS (lon, lat) coordinates and the azimuth angle of the TX
%     antenna's pointing direction in degree, respectively. Note that we
%     have azimuth angles with counter-clockwise as the positive direction
%     starting from +y.
%   - beamwidthInDeg
%     The TX antenna beamwidth in degree.
%   - beamShapeSideLengthInM
%     The length of the two sides originating from the TX location for the
%     polyshape to be generated in meter.
%
% Output:
%   - beamLonLatPolyshape
%     The resulting polyshape with vertices in the form of (lon, lat).
%
% Yaguang Zhang, Purdue, 10/02/2019

% Assume first the TX is located at the UTM system origin and pointing at
% txAzInDeg, find the other two vertices.
v1AngleInDeg = 90+txAzInDeg + beamwidthInDeg/2;
v2AngleInDeg = 90+txAzInDeg - beamwidthInDeg/2;
v1X = cosd(v1AngleInDeg).*beamShapeSideLengthInM;
v2X = cosd(v2AngleInDeg).*beamShapeSideLengthInM;
v1Y = sind(v1AngleInDeg).*beamShapeSideLengthInM;
v2Y = sind(v2AngleInDeg).*beamShapeSideLengthInM;

% Consider the TX UTM location.
[txX, txY, txZone] = deg2utm(txLat, txLon);
v1X = v1X+txX;
v2X = v2X+txX;
v1Y = v1Y+txY;
v2Y = v2Y+txY;

% Convert back to GPS (lat, lon).
[v1Lat, v1Lon] = utm2deg(v1X, v1Y, txZone);
[v2Lat, v2Lon] = utm2deg(v2X, v2Y, txZone);

% Build the polyshape.
beamLonLatPolyshape = polyshape([txLon, txLat; v1Lon, v1Lat; ...
    v2Lon, v2Lat; txLon, txLat]);

end
% EOF