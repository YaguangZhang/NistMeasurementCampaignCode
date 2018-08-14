function [pixelXs, pixelYs, utmZone] = latLon2PixIndices(lats, lons, ...
    LAT_RANGE, LON_RANGE, IMG_RESOLUTION)
%LATLON2PIXINDICES Convert (lat, lon) pairs to the pixel index pairs on the
%vegetation area image.
%
% The centers of the corner pixels will be the boundary points defined by
% element pairs from LAT_RANGE and LON_RANGE. We will also convert GPS
% (lat, lon) pairs to UTM system (x, y) when necessary so that the image
% represents a map in the UTM system.
%
% Inputs:
%   - lats, lons
%     Column vectors. The (lat, lon) pairs to convert.
%   - LAT_RANGE, LON_RANGE
%     A two-element (increasing) column vector to indicate the GPS
%     locations of the boundary pixels.
%   - IMG_RESOLUTION
%     A two-element column vector (width, height) for specifying the
%     resolution of the image.
%
% Outputs:
%   - pixelXs, pixelYs
%     The pixel indices for the input location in the image. Note that we
%     have sticked to the tradition of image processing for indexing pixels
%     (with y indices increasing from top to bottom and the top-left pixel
%     being the origin).
%   - utmZone
%     The UTM zone in which the boundry points lay.
%
% Yaguang Zhang, Purdue, 08/06/2018

numInputPts = length(lats);
[pixelXs, pixelYs] = deal(nan(numInputPts, 1));
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

boundPixXY = nan(4,2);
for idxBoundPt = 1:4
    [curBoundX, curBoundY, curBoundZone] = deg2utm( ...
        boundPixLonLat(idxBoundPt, 2), boundPixLonLat(idxBoundPt, 1));
    assert(strcmp(utmZone, curBoundZone), ...
        'Multiple UTM zones involved for the boundary!');
    
    boundPixXY(idxBoundPt, :) = [curBoundX, curBoundY];
end
minBoundX = boundPixXY(1,1);
maxBoundX = boundPixXY(4,1);
deltaBoundX = maxBoundX - minBoundX;
minBoundY = boundPixXY(1,2);
maxBoundY = boundPixXY(2,2);
deltaBoundY = maxBoundY - minBoundY;

% Close the boundary polygon.
boundPixXYBoundPoly = [boundPixXY; boundPixXY(1,:)];

% Treat input GPS points one by one.
for idxInputPt = 1:numInputPts
    curLat = lats(idxInputPt);
    curLon = lons(idxInputPt);
    [curX, curY, curZone] = deg2utm(curLat, curLon);
    
    % Only return valid outputs if the input point is in the same UTM zone
    % of the boundary and also contained by the boundary polygon.
    if strcmp(utmZone, curZone) && ...
            inpolygon(curX, curY, ...
            boundPixXYBoundPoly(:,1), ...
            boundPixXYBoundPoly(:,2))
        pixelXs(idxInputPt) = ...
            round((imageWidth-1).*(curX-minBoundX)./deltaBoundX+1);
        pixelYs(idxInputPt) = ...
            round((imageHeight-1).*(maxBoundY-curY)./deltaBoundY+1);
    end
end

end

% EOF
