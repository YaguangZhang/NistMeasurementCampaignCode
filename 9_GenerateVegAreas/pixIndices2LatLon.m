function [lats, lons, utmZone] = pixIndices2LatLon(pixelXs, pixelYs, ...
    LAT_RANGE, LON_RANGE, IMG_RESOLUTION)
%PIXINDICES2LATLON Convert the pixel index pairs on the vegetation area
%image to (lat, lon) pairs. This is the convert function of
%latLon2PixIndices.
%
% Please see latLon2PixIndices.m for more information.
%
% Yaguang Zhang, Purdue, 08/06/2018

numInputPairs = length(pixelXs);
[lats, lons] = deal(nan(numInputPairs, 1));
imageWidth = IMG_RESOLUTION(1);
imageHeight = IMG_RESOLUTION(2);

% Get (lat, lon) boundary points. We will start with the top-left piont
% going clock-wise.
boundPixLonLat = [min(LON_RANGE) min(LAT_RANGE); ... % (lon, lat)
    min(LON_RANGE) max(LAT_RANGE); ...
    max(LON_RANGE) max(LAT_RANGE); ...
    max(LON_RANGE) min(LAT_RANGE)];

% Convert all relavent GPS points to UTM coordinates.
[~,~,utmZone] = deg2utm( ...
    boundPixLonLat(1, 2), boundPixLonLat(1, 1));
assert(ischar(utmZone) && (~isempty(utmZone)), ...
    'The UTM zone for boundary is not valid!');

boundXY = nan(4,2);
for idxBoundPt = 1:4
    [curBoundX, curBoundY, curBoundZone] = deg2utm( ...
        boundPixLonLat(idxBoundPt, 2), boundPixLonLat(idxBoundPt, 1));
    assert(strcmp(utmZone, curBoundZone), ...
        'Multiple UTM zones involved for the boundary!');
    
    boundXY(idxBoundPt, :) = [curBoundX, curBoundY];
end
minBoundX = boundXY(1,1);
maxBoundX = boundXY(4,1);
deltaBoundX = maxBoundX - minBoundX;
minBoundY = boundXY(1,2);
maxBoundY = boundXY(2,2);
deltaBoundY = maxBoundY - minBoundY;

% Treat input pixel index pairs one by one.
for idxInputIdxPair = 1:numInputPairs
    curPixX = pixelXs(idxInputIdxPair);
    curPixY = pixelYs(idxInputIdxPair);
    
    % Only return valid outputs if the input index pair is valid and
    % contained by the image.
    if (floor(curPixX)==curPixX) && (floor(curPixX)==curPixX) ...
            && (curPixX>=1) && (curPixY>=1) ...
            && (curPixX<=IMG_RESOLUTION(1)) && (curPixY<=IMG_RESOLUTION(2))
        curX = (curPixX-1)./(imageWidth-1).*deltaBoundX+minBoundX;
        curY = maxBoundY-(curPixY-1)./(imageHeight-1).*deltaBoundY;
        [lats(idxInputIdxPair),lons(idxInputIdxPair)] ...
            = utm2deg(curX,curY,utmZone);
    end
end


end
% EOF