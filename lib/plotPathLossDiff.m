function [hFigDiff] = plotPathLossDiff(pathLossDiffXYZs, dataLabel, ...
    colorRange, rxXY, treeXYs, ...
    flagGenFigSilently)
%PLOTPATHLOSSDIFF Show the path loss difference on a 2D plot with TX and
%trees.
%
% We will show the absolute values of the difference via plot3k and color
% the back ground by the sign (+/-) of the difference via surf.
%
% Inputs:
%   - pathLossDiffXYZs
%     A matrix with each row being (x, y, path loss diff).
%   - colorRange
%     The color range to feed into plot3k.
%   - rxXY
%     The (x, y) coordinates for the RX.
%   - treeXYs
%     A matrix with each row being (x, y) coordinates of a tree to show.
%   - flagGenFigSilently
%     A flag to control whether to show the figure or not.
%
% Yaguang Zhang, Purdue, 10/02/2019

hFigDiff = figure('visible', ~flagGenFigSilently);
hold on; curColormap = colormap('hot'); colormap(curColormap);
plot3k([pathLossDiffXYZs(:, 1:2), abs(pathLossDiffXYZs(:,3))], ...
    'Labels', {'', 'Longitude', 'Latitude', '', ...
    ['abs(', dataLabel, ' - Sim) (dB)']}, ...
    'ColorRange', colorRange);

curZs = pathLossDiffXYZs(:,3);
higherLayerZ = max(abs(curZs))+1;

hTx = plot3(rxXY(1), rxXY(2), higherLayerZ, '^g', 'LineWidth', 2);
[numOfTrees, ~] = size(treeXYs);
hTrees = plot3(treeXYs(:,1), treeXYs(:,2), ...
    ones(numOfTrees).*higherLayerZ, 'b*');
axis tight; view(2);
axisToSet = axis;

%% Background for +/- Sign

% Create meshgrid for surf. We will increase the grid density by a factor
% of 10 to better show the blockage areas.
sufNumPtsPerSide = 1000;
xs = pathLossDiffXYZs(:, 1);
ys = pathLossDiffXYZs(:, 2);
zs = pathLossDiffXYZs(:, 3);

% Set blockage area to 1 and other areas to 0.
boolsIsPos = zs>0;
zs(boolsIsPos) = 1;
zs(~boolsIsPos) = 0;

% Find the ranges for the boundary of interet (BoI) to build a new grid for
% showing the results.
k = convhull([xs, ys]);
xsBoI = xs(k);
ysBoI = ys(k);
xMinBoI = min(xsBoI);
xMaxBoI = max(xsBoI);
yMinBoI = min(ysBoI);
yMaxBoI = max(ysBoI);

[xsNew, ysNew] = meshgrid( ...
    linspace(xMinBoI, xMaxBoI, sufNumPtsPerSide), ...
    linspace(yMinBoI, yMaxBoI, sufNumPtsPerSide));
zsNew = griddata(xs,ys,zs,xsNew,ysNew,'Nearest');

% Ignore points out of the area of interest by seting the z values for them
% to NaN.
[in,on] = inpolygon(xsNew(:), ysNew(:), xsBoI, ysBoI);
boolsPtsToIgnore = ~(in|on);
if any(boolsPtsToIgnore)
    zsNew(boolsPtsToIgnore) = nan;
end

csNew = ones(sufNumPtsPerSide, sufNumPtsPerSide, 3);
csNew = csNew.*repmat(zsNew, [1,1,3]);
% Plot the blockage areas.
surf(xsNew, ysNew, zsNew, csNew, ...
    'FaceAlpha', 0.5, 'EdgeColor', 'none');

hPos = plot(polyshape(nan(3,2)), 'FaceColor', 'w');
hNeg = plot(polyshape(nan(3,2)), 'FaceColor', 'k');

%% Extra Settings
legend([hTx, hTrees(1), hPos, hNeg], ...
    'TX', 'Trees', 'Diff>0', 'Diff<=0', ...
    'Location', 'southeast');
% plotGoogleMapAfterPlot3k(hFigOverviewForSimPLsForGrid, 'satellite');
title(['RMSE = ', num2str(sqrt(mean(curZs(~isnan(curZs)).^2))), ' dB']);
axis(axisToSet);

end
% EOF