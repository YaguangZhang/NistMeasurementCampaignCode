function [foliageDepth, losDist] = estimateFoliageDepthGps( ...
    srcGpsPt, dstGpsPt, vegAreas, VEG_AREA_IMG_META, FLAG_DEBUG)
%ESTIMATEFOLIAGEDEPTHUTM Estimate the folaige depth between two GPS points
%(lat, lon) within the same UTM zone.
%
% Inputs:
%   - srcGpsPt, dstGpsPt
%     The source and destination points (2x1 vectors).
%   - vegAreas, VEG_AREA_IMG_META
%     The image structure reprenting the vegetation area and the meta
%     information for it, e.g. generated by:
%         9_GenerateVegAreas/generateVegAreas.m
%   - FLAG_DEBUG
%     An optional boolean flag. Set this to be true to generate debugging
%     figures.
%
% Output:
%   - foliageDepth
%     The estimated foliage depth in meters.
%
% Yaguang Zhang, Purdue, 08/13/2018

% No need to debug by default.
if ~exist('FLAG_DEBUG', 'var')
    FLAG_DEBUG = false;
end

% Change input GPS points to UTM and verify that they are within the same
% zone.
[inputPtXs, inputPtYs, inputPtZones] = deg2utm( ...
    [srcGpsPt(1); dstGpsPt(1)], [srcGpsPt(2); dstGpsPt(2)]);
assert(strcmp(inputPtZones(1, :), inputPtZones(2, :)), ...
    'Input GPS points should be in the same UTM zone!');

% Calculate the LoS distance.
srcUtmPt = [inputPtXs(1); inputPtYs(1)];
dstUtmPt = [inputPtXs(2); inputPtYs(2)];
losDist = sqrt(sum((srcUtmPt-dstUtmPt).^2));

% Find the pixels corresponding to the input GPS points in the vegetation
% area image.
[srcPixPtX, srcPixPtY] = latLon2PixIndices(srcGpsPt(1), srcGpsPt(2), ...
    VEG_AREA_IMG_META.LAT_RANGE, VEG_AREA_IMG_META.LON_RANGE, ...
    VEG_AREA_IMG_META.IMG_RESOLUTION);
[dstPixPtX, dstPixPtY]  = latLon2PixIndices(dstGpsPt(1), dstGpsPt(2), ...
    VEG_AREA_IMG_META.LAT_RANGE, VEG_AREA_IMG_META.LON_RANGE, ...
    VEG_AREA_IMG_META.IMG_RESOLUTION);

% Find vegetation area pixels along the LoS line segment.
boolsInVegArea = improfile(vegAreas, ...
    [srcPixPtX, dstPixPtX], [srcPixPtY, dstPixPtY])>0;

% Calculate the in-vegetation-area ratio.
inVegAreaRatio = sum(boolsInVegArea)/length(boolsInVegArea);

% Estimate the LoS foliage depth accordingly.
foliageDepth = losDist.*inVegAreaRatio;

% For debugging.
if FLAG_DEBUG
    % Source and destination locations on a map.
    figure; hold on;
    hSrc = plot(srcGpsPt(2), srcGpsPt(1), '^g');
    hDst = plot(dstGpsPt(2), dstGpsPt(1), '*r');
    plot([srcGpsPt(2);dstGpsPt(2)], [srcGpsPt(1);dstGpsPt(1)], 'b-');
    plot_google_map('MapType', 'satellite');
    legend([hSrc,hDst], 'Source', 'Destination');
    
    % Vegetation image with the source and destination locations.
    figure; hold on;
    [pixelIndicesX, pixelIndicesY] ...
        = meshgrid(1:VEG_AREA_IMG_META.IMG_RESOLUTION(1), ...
        -1:-1:-VEG_AREA_IMG_META.IMG_RESOLUTION(2));
    boolsIsVegPixel = vegAreas(:)>0;
    plot(pixelIndicesX(boolsIsVegPixel), ...
        pixelIndicesY(boolsIsVegPixel), '.k');
    hSrc = plot(srcPixPtX, -srcPixPtY, '^g');
    hDst = plot(dstPixPtX, -dstPixPtY, '*r');
    plot([srcPixPtX;dstPixPtX], [-srcPixPtY;-dstPixPtY], 'b-');
    legend([hSrc,hDst], 'Source', 'Destination');
end

end
% EOF