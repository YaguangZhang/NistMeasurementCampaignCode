function [ firstFresZoneArea ] ....
    = exploreForPixelsInFirstFres( ...
    firstFresZoneArea, startPixelCoors, ...
    srcUtmPt3D, dstUtmPt3D, ...
    VEG_AREA_IMG_META, fCarrierGHz)
%EXPLOREFORPIXELSINFIRSTFRES Explore up, down, left, and right for pixels
%in the 1s Fresnel zone for the foliage analysis.
%
% We will explore the four directions seperately. For each direction, we
% will go as far as we can until a pixel that is out of the zone is
% encountered.
%
% Inputs:
%   - firstFresZoneArea
%     The known mapping for the pixels. 1: It is a pixel in the 1st Fresnel
%     zone.
%   - startPixelCoors
%     The pixel coordinates for the exploration starting point, in the form
%     of [cx, cy]. Note that we have followed the image processing
%     convention for indexing.
%   - srcUtmPt3D, dstUtmPt3D
%     The source and destination points (3x1 vectors) in (x, y, alt).
%   - VEG_AREA_IMG_META
%     The meta information for the image structure reprenting the
%     vegetation area, e.g. generated by:
%         9_GenerateVegAreas/generateVegAreas.m
%   - fCarrierGHz
%     The signal frequency in GHz.
%
% Output:
%   - firstFresZoneArea
%     The updated boolean matrix indicating which pixels are in the 1st
%     Fresnel zone.
%
% Yaguang Zhang, Purdue, 09/24/2018

% Explore up.
[ firstFresZoneArea ] ....
    = exploreForPixelsInFirstFresOneDir( ...
    firstFresZoneArea, startPixelCoors, [0,-1], srcUtmPt3D, dstUtmPt3D, ...
    VEG_AREA_IMG_META, fCarrierGHz);
% Explore left.
[ firstFresZoneArea ] ....
    = exploreForPixelsInFirstFresOneDir( ...
    firstFresZoneArea, startPixelCoors, [-1,0], srcUtmPt3D, dstUtmPt3D, ...
    VEG_AREA_IMG_META, fCarrierGHz);
% Explore down.
[ firstFresZoneArea ] ....
    = exploreForPixelsInFirstFresOneDir( ...
    firstFresZoneArea, startPixelCoors, [0,1], srcUtmPt3D, dstUtmPt3D, ...
    VEG_AREA_IMG_META, fCarrierGHz);
% Explore right.
[ firstFresZoneArea ] ....
    = exploreForPixelsInFirstFresOneDir( ...
    firstFresZoneArea, startPixelCoors, [1,0], srcUtmPt3D, dstUtmPt3D, ...
    VEG_AREA_IMG_META, fCarrierGHz);

end

% EOF