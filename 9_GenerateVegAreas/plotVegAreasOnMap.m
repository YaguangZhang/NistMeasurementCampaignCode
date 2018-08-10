function [ hVegAreas ] = plotVegAreasOnMap( ...
    hFig, vegAreas, VEG_AREA_IMG_META, vegAreaZ)
%PLOTVEGAREASONMAP Plot the vegetation pixels as transparent dots on the
%map specified by hInterVegAreaMarker.
%
% Inputs:
%   - hFig
%     The handle to the figure to add the plot.
%   - vegAreas, VEG_AREA_IMG_META
%     The image to represent the vegetation area and the corresponding meta
%     data for it, generated by generateVegAreas.mat.
%   - vegAreaZ
%     The z value that will be used to plot the vegetation area.
%
% Output:
%   - hVegAreas
%     The handle to the vegetation area plot.
%
% Yaguang Zhang, Purdue, 08/08/2018

alphaForTransparency = 0.3;
imageResolution = size(vegAreas);

figure(hFig); hold on;

hVegAreas = surf( ...
    VEG_AREA_IMG_META.LONS, ...
    VEG_AREA_IMG_META.LATS, ...
    ones(imageResolution).*vegAreaZ, ...
    'FaceAlpha','flat',...
    'AlphaData', vegAreas.*alphaForTransparency, ...
    'AlphaDataMapping', 'none', ...
    'FaceColor','w','LineStyle','none');

end

% EOF