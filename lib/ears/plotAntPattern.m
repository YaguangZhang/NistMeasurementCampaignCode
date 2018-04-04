function [ hPat2DRef, hPat3DRef, ...
    hInterPat3DOnLineLinear, hInterPat3DOnLineDb, ...
    hInterPat3DWeightedSumLinear, hInterPat3DWeightedSumDb ] ...
    = plotAntPattern( patAz, patEl, absPathWithPrefixToSavePlots )
%PLOTANTPATTERN Plot the antenna pattern specified by the inputs patAz and
%patEl.
%
% Inputs:
%   - patAz, patEl
%     The antenna patterns, for the Azimuth and Elevation sweeps,
%     respectively; Each of which is a struct containing fields:
%       - azs
%         The azimuth angles in degree from set [0, 360).
%       - els
%         The elevation angles in degree from set [0, 360).
%       - amps
%         The linear amplitudes of the samples.
%       - phases
%         The phases of the samples.
%     All of these fields contains a column vector with each row
%     corresponding to a sweep sample.
%   - absPathWithPrefixToSavePlots
%     A optional string to specify where to save the plots. When it is not
%     present, the plots will only be generated without saving.
%
% For the reference input sweep data, both a 2D illsutration and a 3D
% illustrations will be generated. For the interpolated data, a 3D
% illustration will be generated for each interpolation method.
%
% Ref:
% http://antennatutorials.blogspot.com/2013/05/radiation-pattern-of-half-wave-dipole.html
%
% Yaguang Zhang, Purdue, 10/02/2017

%% 2D Plot
% Update: We will also plot the HPBW for both cases on the angular plane.

hPat2DRef = figure('units','normalized', ...
    'outerposition',[0.1 0.05 0.8 0.9], 'Name','hPat2DRef');
% Azimuth.
polarAx = subplot(2,2,1,polaraxes);
[~, hpbw] = plotAntPlanePat(polarAx, patAz.azs, antPatLinearToDb(patAz.amps));
title({'Azimuth Plane Pattern (Relative to the Minimum Amplitude)'; ...
    ['HPBW = ', num2str(hpbw), ' degrees']});

subplot(2,2,3);
angles = patAz.azs;
plot(angles, antPatLinearToDb(patAz.amps));
title('Azimuth Sweep Data'); grid minor;
curAxis = axis; axis([min(angles), max(angles), curAxis(3:4)]);
xlabel('Azimuth'); ylabel('Normalized Amplitude (dB)');

% Elevation.
polarAx = subplot(2,2,2,polaraxes);
[~, hpbw] = plotAntPlanePat(polarAx, patEl.els, antPatLinearToDb(patEl.amps));
title({'Elevation Plane Pattern (Relative to the Minimum Amplitude)'; ...
    ['HPBW = ', num2str(hpbw), ' degrees']});

subplot(2,2,4);
angles = patEl.els;
plot(angles, antPatLinearToDb(patEl.amps));
title('Elevation Sweep Data'); grid minor;
curAxis = axis; axis([min(angles), max(angles), curAxis(3:4)]);
xlabel('Elevation'); ylabel('Normalized Amplitude (dB)');

% Also generate dedicated plots for the antenna patterns.
[~, ~, hPat2DAz] = plotAntPlanePat('', patAz.azs, ...
    antPatLinearToDb(patAz.amps), true);
[~, ~, hPat2DEl] = plotAntPlanePat('', patEl.els, ...
    antPatLinearToDb(patEl.amps), true);

%% 3D Plot

% Plot the interpolation results.
numPtsPerDim = 1000;
hInterPat3DOnLineLinear = plotInterPat3D( patAz, patEl, ...
    'OnLine', false, numPtsPerDim);
set(hInterPat3DOnLineLinear, 'Name', 'interPat3DOnLineLinear');
hInterPat3DOnLineDb = plotInterPat3D( patAz, patEl, ...
    'OnLine', true, numPtsPerDim);
set(hInterPat3DOnLineDb, 'Name', 'interPat3DOnLineDb');
hInterPat3DWeightedSumLinear = plotInterPat3D( patAz, patEl, ...
    'WeightedSum', false, numPtsPerDim);
set(hInterPat3DWeightedSumLinear, 'Name', 'interPat3DWeightedSumLinear');
hInterPat3DWeightedSumDb = plotInterPat3D( patAz, patEl, ...
    'WeightedSum', true, numPtsPerDim);
set(hInterPat3DWeightedSumDb, 'Name', 'interPat3DWeightedSumDb');

% We will also plot the sweep data, just for reference.
hPat3DRef = plotRefPat3D( patAz, patEl);

%% Save the Plots
if nargin>2
    absPathCurFile = [absPathWithPrefixToSavePlots, 'pat2DRef'];
    saveas(hPat2DRef, [absPathCurFile, '.png']);
    saveas(hPat2DRef, [absPathCurFile, '.fig']);
    
    absPathCurFile = [absPathWithPrefixToSavePlots, 'pat3DRef'];
    saveas(hPat3DRef, [absPathCurFile, '.png']);
    saveas(hPat3DRef, [absPathCurFile, '.fig']);
    
    absPathCurFile = [absPathWithPrefixToSavePlots, ...
        'interPat3DOnLineLinear'];
    saveas(hInterPat3DOnLineLinear, [absPathCurFile, '.png']);
    saveas(hInterPat3DOnLineLinear, [absPathCurFile, '.fig']);
    
    absPathCurFile = [absPathWithPrefixToSavePlots, 'interPat3DOnLineDb'];
    saveas(hInterPat3DOnLineDb, [absPathCurFile, '.png']);
    saveas(hInterPat3DOnLineDb, [absPathCurFile, '.fig']);
    
    absPathCurFile = [absPathWithPrefixToSavePlots, ...
        'interPat3DWeightedSumLinear'];
    saveas(hInterPat3DWeightedSumLinear, [absPathCurFile, '.png']);
    saveas(hInterPat3DWeightedSumLinear, [absPathCurFile, '.fig']);
    
    absPathCurFile = [absPathWithPrefixToSavePlots, 'interPat3DWeightedSumDb'];
    saveas(hInterPat3DWeightedSumDb, [absPathCurFile, '.png']);
    saveas(hInterPat3DWeightedSumDb, [absPathCurFile, '.fig']);
    
    absPathCurFile = [absPathWithPrefixToSavePlots, 'hPat2DAz'];
    saveas(hPat2DAz, [absPathCurFile, '.png']);
    saveas(hPat2DAz, [absPathCurFile, '.fig']);
    
    absPathCurFile = [absPathWithPrefixToSavePlots, 'hPat2DEl'];
    saveas(hPat2DEl, [absPathCurFile, '.png']);
    saveas(hPat2DEl, [absPathCurFile, '.fig']);
end

end
% EOF