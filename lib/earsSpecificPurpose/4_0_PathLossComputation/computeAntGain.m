function [ gain, hDebugFig, hDebugFigMap] ...
    = computeAntGain(lat0, lon0, alt0, ...
    az0, el0, ...
    lat, lon, alt, ...
    antPatAz, antPatEl, FLAG_DEBUG)
%COMPUTEANTGAIN Compute the antenna gain for the antenna located at (lat0,
%lon0, alt0) facing at (az0, el0), with its pattern specified by the
%Azimuth and Elevation sweep results (antPatAz, antPatEl), from a
%transmitter / to a receiver located at (lat, lon, alt).
%
%   - antPatAz, antPatEl
%     The antenna patterns, for the Azimuth and the Elevation sweeps; Each
%     of them is a struct containing fields:
%       - azs
%         The azimuth angles in degree.
%       - els
%         The elevation angles in degree.
%       - amps
%         The linear amplitudes of the samples.
%       - phases
%         The phases of the samples.
%
% Yaguang Zhang, Purdue, 10/04/2017

% Set this to be true to show a fully interpolated 3D antenna pattern;
% Otherwise, only the reference azimuth and elevation patterns will be
% shown.
FLAG_3D_ATN_PAT = true;

if nargin<11
    FLAG_DEBUG = false;
end

hDebugFig = nan;
hDebugFigMap = nan;

% Set this to true to also show all the intermediant vectors for debugging.
FLAG_VERBOSE = false;

% Convert GPS to UTM so that the location coordinates will all in the unit
% of meter.
[x0Easting, y0Northing, zone0] = deg2utm(lat0, lon0);
[xEasting, yNorthing, zone] = deg2utm(lat, lon);

% Note that from deg2utm we have +x as easting and +y as northing (for
% Azimuth), but Matlab computes Azimuth from +x.
x0 = y0Northing;
y0 = -x0Easting;
x = yNorthing;
y = -xEasting;

if ~strcmp(zone0, zone)
    error('computeAntGain: The two locations specified are not in the same UTM zone!');
end

az0R = deg2rad(az0);
el0R = deg2rad(el0);

% Compute the direction vector from (x0,y0,alt0) to (x,y,alt).
xDiff = x-x0;
yDiff = y-y0;
zDiff = alt-alt0;
% Compute the (az,el,r) for that direction. Note that the resulted angles
% are in radians.
[az1R, el1R, r1] ... % Unit: radian.
    = cart2sph(xDiff, yDiff, zDiff);

% Get the spherical coordinates for that direction with respect of the
% origin antenna's point of view. One can think of the operation as:
%   1. Plot the vector for the direction from the origin antenna to the
%   destination antenna, which is (xDiff, yDiff, zDiff) here;
%  2.  Plot in the same Cartesian coordinate system a vector for the
%  direction to which origin antenna is facing; We essentially use
%  sph2cart(az0, el0, 1) here.
%   3. Rotate these two vectors at the same time, first along Azimuth then
%   along Elevation (of vector in (2)), until the vector in (2) align
%   perfectly with (1,0,0), i.e. the position x axis.
%  4.  The resulted vector in (1) will give us the (az, el) needed.

% Rotate along Azimuth.
[x1, y1, z1] = sph2cart(az1R-deg2rad(az0), el1R, r1);

% Rotate along Elevation of vector in (2). Essentially, we keep y1
% untouched, and rotate (x1, 0, z1) to get the new x1 and z1.

%This method won't always work.
% [az2R, el2R, r2] ... % Unit: radian.
%     = cart2sph(x1, 0, z1);
% [x1F, ~, z1F] = sph2cart(az2R, el2R-deg2rad(el0), r2);

% We will use rotation matrix instead. Note we need the transformation
% matrix to convert a vector in our normal (x, y, z) coordinate system to
% the rotated system.
rotMatR = nan(3,3);
[rotMatR(1,1), rotMatR(1,2), rotMatR(1,3)] = sph2cart(0, -el0R, 1);
rotMatR(2,:) = [0 1 0];
[rotMatR(3,1), rotMatR(3,2), rotMatR(3,3)] = sph2cart(0, -el0R + pi/2, 1);
rotResult = [x1, y1, z1] * rotMatR;
x1F = rotResult(1);
z1F = rotResult(3);

% Get the final result.
[azR, elR, ~] = cart2sph(x1F, y1, z1F);

% Convert them to degree as required by antPatInter.m.
az = rad2deg(azR);
el = rad2deg(elR);
% Now get the antenna gain from the interpolated patten. We will
% interpolate directly using the linear amplitudes of the antenna pattern
% with the "WeightedSum" method.
if FLAG_DEBUG
    FLAG_INTER_IN_DB = false;
    numPtsPerDim = 1000; % Controls the interpolation density.
end
INTER_METHOD = 'WeightedSum';
amps = antPatInter(antPatAz, antPatEl, az, el, 'WeightedSum');

gain = antPatLinearToDb(amps);

% Plot the positions.
if FLAG_DEBUG
    % For scaling the plotted components.
    gridStep = 1; % In meter.
    gridWidth = sqrt(xDiff^2+yDiff^2)/2; % In meter.
    
    % The origin antenna orientation.
    [xOri, yOri, zOri] = sph2cart(az0R, el0R, gridWidth/2);
    newBasisX = [xOri; yOri; zOri];
    newBasisX = newBasisX/norm(newBasisX);
    xNewX = x0+xOri;
    yNewX = y0+yOri;
    zNewX = alt0+zOri;
    % The other two axes for the new coordinate system.
    [xNewYDiff, yNewYDiff, zNewYDiff] ...
        = sph2cart(az0R+pi/2, 0, gridWidth/2);
    newBasisY = [xNewYDiff; yNewYDiff; zNewYDiff];
    newBasisY = newBasisY/norm(newBasisY);
    xNewY = x0+xNewYDiff;
    yNewY = y0+yNewYDiff;
    zNewY = alt0+zNewYDiff;
    [xNewZDiff, yNewZDiff, zNewZDiff] = sph2cart( ...
        az0R+pi, pi/2-el0R, gridWidth/2);
    newBasisZ = [xNewZDiff; yNewZDiff; zNewZDiff];
    newBasisZ = newBasisZ/norm(newBasisZ);
    xNewZ = x0+xNewZDiff;
    yNewZ = y0+yNewZDiff;
    zNewZ = alt0+zNewZDiff;
    
    % The new X-Y plane.
    p1 = [x0, y0, alt0];
    p2 = [xNewY, yNewY, zNewY];
    p3 = [xNewX, yNewX, zNewX];
    normal = cross(p1 - p2, p1 - p3);
    d = p1(1)*normal(1) + p1(2)*normal(2) + p1(3)*normal(3);
    d = -d;
    xGrid = (-(gridWidth/2):gridStep:(gridWidth/2)) + x0;
    yGrid = (-(gridWidth/2):gridStep:(gridWidth/2)) + y0;
    [X,Y] = meshgrid(xGrid,yGrid);
    Z = (-d - (normal(1)*X) - (normal(2)*Y))/normal(3);
    
    % Antenna pattern.
    if FLAG_3D_ATN_PAT
        [ antPatPtsXs,antPatPtsYx,antPatPtsZs, ...
            antPatAmpInDb, antPatAmpInDbRel ] ...
            = interpAntPatsIn3D( ...
            antPatAz, antPatEl, ...
            numPtsPerDim, ...
            FLAG_INTER_IN_DB, INTER_METHOD, true);
        % Rearrange the outputs to vector columns and scale the coordinates
        % to a max length.
        antPatAmpInDb = antPatAmpInDb(1:end)';
        [ maxAntPatAmpInDbRel ] = max(antPatAmpInDbRel(1:end));
        scalingFactor = (gridWidth/2*0.8)/maxAntPatAmpInDbRel;
        antPatPts = [antPatPtsXs(1:end); ...
            antPatPtsYx(1:end); ...
            antPatPtsZs(1:end)]' ...
            .*scalingFactor;
    else
        antPatAmpInDb = antPatLinearToDb([antPatAz.amps; antPatEl.amps]);
        antPatRelAmpInDb = antPatAmpInDb-min(antPatAmpInDb);
        [antPatPtsXs, antPatPtsYx, antPatPtsZs] ...
            = sph2cart(deg2rad([antPatAz.azs; antPatEl.azs]), ...
            deg2rad([antPatAz.els; antPatEl.els]), antPatRelAmpInDb);
        antPatPts = [antPatPtsXs, antPatPtsYx, antPatPtsZs];
    end
    
    % Convert the pattern to the new coordinate system.
    newCoorMatR = [newBasisX, newBasisY, newBasisZ]';
    antPatPtsNewCoor = bsxfun(@plus, antPatPts*newCoorMatR, ...
        [x0, y0, alt0]);
    
    hDebugFig = figure('units','normalized', ...
        'outerposition',[0.1 0.05 0.8 0.9]);
    % The UTM illustration with all the details.
    hold on; colormap jet;
    % The origin and destination.
    hOri = plot3(x0, y0, alt0, '^w', 'MarkerFaceColor', 'red');
    hDes = plot3(x, y, alt, 'sw', 'MarkerFaceColor', 'blue');
    % Link them with a gray line.
    plot3([x0; x], [y0, y], [alt0, alt], '--', ...
        'Color', ones(1,3)*0.7);
    % Plot the origin antenna orientation.
    axisLineWidth = 2;
    hOriOri = plot3([x0, xNewX], [y0, yNewX], [alt0, zNewX], ...
        '-r', 'LineWidth', axisLineWidth);
    % Plot the Cartesian coordinate system defined by the origin antenna
    % orientation.
    hNewY = plot3([x0, xNewY], [y0, yNewY], ...
        [alt0, zNewY], '-g', 'LineWidth', axisLineWidth);
    hNewZ = plot3([x0, xNewZ], [y0, yNewZ], ...
        [alt0, zNewZ], '-b', 'LineWidth', axisLineWidth);
    % Plot the new x-y plane. We use plot3 instead of mesh because the
    % later one will mess up the color bar for plot3k.
    hNewXY = plot3(X(:), Y(:), Z(:), '.', 'Color', ones(1,3)*0.9);
    % Plot the reference antenna pattern.
    plot3k(antPatPtsNewCoor, 'ColorData', antPatAmpInDb, ...
        'ColorRange', [min(antPatAmpInDb) max(antPatAmpInDb)]);
    if FLAG_VERBOSE
        plot3(x0+x1, y0+y1, alt0+z1, 'k+');
        plot3(x0+x1F, y0+y1, alt0+z1F, 'k*');
    end
    % Embed the final results (az, el) in the figure title for human
    % inspectation.
    title({['Origin (with Facing Direction az=',num2str(az0), ...
        ',el=',num2str(el0), ') and Destination']; ...
        ['Computed (az, el) from the Origin''s view: (', ...
        num2str(az), ',', num2str(el),') degrees'];
        ['Gain = ', num2str(gain), ...
        ' dB according to the Colored Antenna Pattern (in dB)']})
    hold off; grid minor;
    view(-105, 20); axis equal;
    xlabel('x (northing in m)'); ylabel('y (westing in m)');
    zlabel('z (altitude in m)');
    legend([hOri, hDes, hOriOri, hNewY, hNewZ, hNewXY], ...
        'Origin', 'Dest', ...
        'Origin Orientation (X)', 'Y', 'Z', 'Plane X-Y');
    
    % The (lat, lon) plot to show the origin and destination on a map.
    hDebugFigMap = figure('units','normalized', ...
        'outerposition',[0.2 0.05 0.6 0.9]);
    hold on;
    hOri = plot3(lon0, lat0, alt0, '^w', 'MarkerFaceColor', 'red');
    hDes = plot3(lon, lat, alt, 'sw', 'MarkerFaceColor', 'blue');
    plot3([lon0; lon], [lat0, lat], [alt0, alt], '--', ...
        'Color', ones(1,3)*0.7);
    hold off; plot_google_map('MapType', 'satellite');
    view(-10, 45); grid minor;
    xlabel('Lon'); ylabel('Lat'); zlabel('Alt (m)');
    legend([hOri, hDes], 'Origin', 'Dest');
end

end

%EOF