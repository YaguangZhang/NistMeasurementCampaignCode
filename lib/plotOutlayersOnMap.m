function hFigOutlayersOnMap = plotOutlayersOnMap( ...
    x, lons, lats, numOfSigmasOutlayer, figAxisToSet)
%PLOTOUTLAYERSONMAP Plot the outlayers of x on a map.
%
% Inputs:
%   - x, lons, lats
%     The data to inspect as well as their GPS locations, respectively.
%   - numOfSigmasOutlayer
%     The number of sigmas to use as the shreshold to determine outlayers.
%     Data points that are out of this many sigma range around the mean
%     will be considered as outlayers.
%   - figAxisToSet
%     An optional vector for setting the axis ranges.
%
% Output:
%   - hFigOutlayersOnMap
%     The handle to the output figure.
%
% Yaguang Zhang, Purdue, 10/17/2019

% By default, find the flag flagGenFigSilently in the base workspace. If
% not found, flagGenFigSilently will be set to false.
try
    flagGenFigSilently = evalin('base', 'flagGenFigSilently');
catch
    flagGenFigSilently = false;
end

if isempty(x)
    hFigOutlayersOnMap = figure('visible', ~flagGenFigSilently);
else
    xSignma = std(x);
    xMean = mean(x);
    
    boolsIsOutlayer = abs(x-xMean)>(numOfSigmasOutlayer.*xSignma);
    
    floatFomatter = '%.2f';
    zShift = -min(x);
    
    hFigOutlayersOnMap = figure('visible', ~flagGenFigSilently); hold on;
    plot3(lons, lats, zeros(size(lons)), '.', 'Color', ones(1,3).*0.7);
    colormap hot;
    [~, ~, hCb] = plot3k([lons(boolsIsOutlayer), lats(boolsIsOutlayer), ...
        x(boolsIsOutlayer)+zShift]);
    hCb.TickLabels = arrayfun(@(idx) ...
        {num2str(str2double(hCb.TickLabels{idx})-zShift)}, ...
        1:length(hCb.Ticks))';
    if exist('figAxisToSet', 'var')
        axis(figAxisToSet);
    end
    view(2); plotGoogleMapAfterPlot3k(hFigOutlayersOnMap, 'satellite');
    xticks([]); yticks([]); xlabel('Longitude'); ylabel('Latitude');
    title(['Outlayers on Map (mu = ',  ...
        num2str(xMean, floatFomatter), ', sigma =', ...
        num2str(xSignma, floatFomatter), ')']);
end
end
% EOF