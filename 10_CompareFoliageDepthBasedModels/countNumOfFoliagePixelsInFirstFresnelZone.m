function [ numOfFoliagePixelsInFirstFresnel, curFirstFresZoneArea] ...
    = countNumOfFoliagePixelsInFirstFresnelZone(...
    srcGpsPt3D, dstGpsPt3D, ...
    vegAreas, VEG_AREA_IMG_META, fCarrierGHz, FLAG_DEBUG)
%COUNTNUMOFFOLIAGEPIXELSINFIRSTFRESNELZONE Count the number of foliage
%pixels in the first Fresnel zone between the source and the destination.
%
% Inputs:
%   - srcGpsPt3D, dstGpsPt
%     The source and destination points (3x1 vectors) in (lat, lon, alt).
%   - vegAreas, VEG_AREA_IMG_META
%     The image structure reprenting the vegetation area and the meta
%     information for it, e.g. generated by:
%         9_GenerateVegAreas/generateVegAreas.m
%   - fCarrierGHz
%     The signal frequency in GHz.
%   - FLAG_DEBUG
%     An optional boolean flag. Set this to be true to generate debugging
%     figures.
%
% Output:
%   - numOfFoliagePixelsInFirstFresnel
%     The resulting number of foliage pixels in the first Fresnel zone
%     between the source and the destination.
%
% Yaguang Zhang, Purdue, 09/21/2018

% We will reuse countNumOfTreesInFirstFresnelZone.
if ~exist('countNumOfTreesInFirstFresnelZone', 'file')
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', ...
        '4_FoliageAttenuationEstimation'));
end

% No need to debug by default.
if ~exist('FLAG_DEBUG', 'var')
    FLAG_DEBUG = false;
end

% Change input GPS points to UTM and verify that they are within the same
% zone.
[inputPtXs, inputPtYs, inputPtZones] = deg2utm( ...
    [srcGpsPt3D(1); dstGpsPt3D(1)], [srcGpsPt3D(2); dstGpsPt3D(2)]);
assert(strcmp(inputPtZones(1, :), inputPtZones(2, :)), ...
    'Input GPS points should be in the same UTM zone!');

% Source and destination locations.
srcUtmPt3D = [inputPtXs(1), inputPtYs(1), srcGpsPt3D(3)];
dstUtmPt3D = [inputPtXs(2), inputPtYs(2), dstGpsPt3D(3)];

% Find the pixels corresponding to the input GPS points in the vegetation
% area image.
[srcPixPtX, srcPixPtY] = latLon2PixIndices( ...
    srcGpsPt3D(1), srcGpsPt3D(2), ...
    VEG_AREA_IMG_META.LAT_RANGE, VEG_AREA_IMG_META.LON_RANGE, ...
    VEG_AREA_IMG_META.IMG_RESOLUTION);
[dstPixPtX, dstPixPtY]  = latLon2PixIndices( ...
    dstGpsPt3D(1), dstGpsPt3D(2), ...
    VEG_AREA_IMG_META.LAT_RANGE, VEG_AREA_IMG_META.LON_RANGE, ...
    VEG_AREA_IMG_META.IMG_RESOLUTION);

% We will build a matrix corresponding to the foliage image showing which
% pixels are in the first Fresnel zone or not.
curFirstFresZoneArea = nan(size(vegAreas));
% Start with pixels on the LoS path.
[cx, cy, ~] = improfile(vegAreas, ...
    [srcPixPtX, dstPixPtX], [srcPixPtY, dstPixPtY]);
% Round the indices to integers.
cx = round(cx); cy = round(cy); numStartPts = length(cx);
% Initialization.
for idxStartPt = 1:numStartPts 
    curFirstFresZoneArea(cy(idxStartPt), cx(idxStartPt)) = 1;
end
% Explore from these starting poinst.
for idxStartPt = 1:numStartPts
    curFirstFresZoneArea ...
        = exploreForPixelsInFirstFres( ...
        curFirstFresZoneArea, [cx(idxStartPt), cy(idxStartPt)], ...
        srcUtmPt3D, dstUtmPt3D, ...
        VEG_AREA_IMG_META, fCarrierGHz);
end
% Unexplored areas are not in the zone.
curFirstFresZoneArea(isnan(curFirstFresZoneArea(:))) = 0;

% Final result.
numOfFoliagePixelsInFirstFresnel = sum( ...
    curFirstFresZoneArea(:) & vegAreas(:));

% Plot for debugging.
if FLAG_DEBUG
    
    if ~exist('plotVegAreasOnMap', 'file')
        addpath(fullfile(fileparts(mfilename('fullpath')), '..', ...
            '9_GenerateVegAreas'));
    end
    
    interceptionAreas = and(curFirstFresZoneArea, vegAreas);
    
    hVegAreaOnMap = figure; hold on;
    hTx = plot3(srcGpsPt3D(2), srcGpsPt3D(1), srcGpsPt3D(3), '^r');
    hRx = plot3(dstGpsPt3D(2), dstGpsPt3D(1), dstGpsPt3D(3), 'xr');
    plot_google_map('MapType' , 'satellite');
    % Color the vegetation areas on the figure.
    hVegAreas = plotVegAreasOnMap(hVegAreaOnMap, vegAreas, ...
        VEG_AREA_IMG_META, 1, 'g');
    % Color the in-zone part.
    hVegAreasInZone = plotVegAreasOnMap(hVegAreaOnMap, ...
        curFirstFresZoneArea, ...
        VEG_AREA_IMG_META, 2, 'y');
    % Color the intersection part.
    hIntersection = plotVegAreasOnMap(hVegAreaOnMap, interceptionAreas, ...
        VEG_AREA_IMG_META, 3, 'r');
    legend([hTx, hRx, hVegAreas, hVegAreasInZone, hIntersection], ...
        'Tx', 'Rx', 'Foliage', '1st Fresnel Zone', 'Intersection', ...
        'Location', 'northwest');
    xlabel('x'); ylabel('y'); zlabel('z');
end

end

% EOF