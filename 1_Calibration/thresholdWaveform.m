function [ waveformThresholded, boolsEliminatedPts, hDebugFig ] ...
    = thresholdWaveform( waveform, flagDebug )
%THRESHOLDWAVEFORM Thresholded the waveform (a real vector) to eliminate
%corss-correlation and system noise.
%
% We have implemented the algorithm provided by Professor Chris Anderson:
%
%   Find the standard deviation in the first 5-10% of samples before the
%   first peak, set the threshold at 3*sigma above 0, and remove everything
%   below the threshold.
%
% For complex signals, you may want to apply this function twice to its
% real and imaginary parts, separately.
%
% Optionally, if flagDebug is set and its value is true, we will generate
% plots for debugging.
%
% Update: we will use 3.28*sigma shifted by the mean to determine the
% threshold, so that more noise will be eliminated.
%
% Yaguang Zhang, Purdue, 08/14/2017

%% Parameters

% Set this to be true to first convert the wave form to dB.
FLAG_PROCESS_IN_DB = true;
% Different ways to find the first present signal / peak.
%   - 'TallEnoughPositive'
%     The first one taller than 90% of the tallest peak.
%   - 'TallEnoughAbs'
%     The first one taller than 90% of the tallest peak in terms of
%     amplitude.
%   - 'TallEnoughAboveMedian'
%     The first one with amplitude taller than 50% of the tallest peak
%     after getting rid of 2 times median of the data.
%   - 'FindPeaks'
%     Built-in Matlab findpeaks function.
PEAK_SEARCH_METHOD = 'TallEnoughAbs';
% For setting the threshold, how many sigmas are used.
try
    NUM_SIGMA_FOR_THRESHOLD = evalin('base', 'NUM_SIGMA_FOR_THRESHOLD');
catch
    warning('    NUM_SIGMA_FOR_THRESHOLD not set in the base workspace. Will use the default value 3 for it.')
    NUM_SIGMA_FOR_THRESHOLD = 3;
end
% Number of samples to discard at the beginning. Note this is a local
% variable and is only used in this function.
numStartSampsToDiscard = 0;
% Maximum number of samples to consider for noise elimination.
maxNumSampsToConsider = 5*10^5;
% Find the first peak after the sample index specified here to make sure
% the noise sigma is computed from enough samples.
idxMinPeak = 2000;
% Relative range for computing noise sigma, from the start of considered
% sample segment to the peak. Ideally, all the samples in the range should
% be noise samples. Here we choose to use the samples from 45% to 55%
% because the part right before the first peak found is more likely to be
% noise.
relativeRangeForNoise = [0.875, 0.975];
% Sample rate used for GnuRadio.
try
    Fs = evalin('base', 'Fs');
catch
    warning('GnuRadio sample frequency Fs not found in the base workspace.')
    warning('Will use the default value 1.04 * 10^6.')
    Fs = 1.04 * 10^6;
end

hDebugFig = nan;
if nargin < 2
    flagDebug = false;
end

%% Algorithm

% Discard the first numFirstSampsToDiscard samples, as well as the samples
% out of the maximum-number-of-samples-to-consider range.
waveformToCons = waveform((numStartSampsToDiscard+1):end);
if length(waveformToCons)>maxNumSampsToConsider
    waveformToCons = waveformToCons(1:maxNumSampsToConsider);
end

waveformToConsAmps = abs(waveformToCons);
% Shift twice median with mean as the bottom line of all the peaks.
twiceMedianShifted = 2*median(waveformToConsAmps)+mean(waveformToConsAmps);

switch PEAK_SEARCH_METHOD
    case 'TallEnoughPositive'
        RATIO_VS_TALLEST = 0.9;
        [ampMax, ~] = max(waveformToCons(idxMinPeak:end));
        idxRefSample = find(waveformToCons(idxMinPeak:end) ...
            >RATIO_VS_TALLEST.*ampMax, 1)+idxMinPeak-1;
        numSubFigs = 2;
    case 'TallEnoughAbs'
        RATIO_VS_TALLEST = 0.9;
        [ampMax, ~] = max(waveformToConsAmps(idxMinPeak:end));
        idxRefSample = find(waveformToConsAmps(idxMinPeak:end) ...
            >RATIO_VS_TALLEST.*ampMax, 1)+idxMinPeak-1;
        numSubFigs = 2;
    case 'TallEnoughAboveMedian'
        % Find the peak that is "tall enough":
        %   1. Use amplitude of the waveform
        %    2. For all amplitudes, minus 2*median
        %   3. Eliminate negative points
        %    4. Find the first sample higher than RATIO_VS_TALLEST of the
        %    tallest
        %   one.
        RATIO_VS_TALLEST = 0.5;
        
        shiftedWaveformAmps = waveformToConsAmps - twiceMedianShifted;
        shiftedWaveformAmps(shiftedWaveformAmps<0) = 0;
        [ampMax, ~] = max(shiftedWaveformAmps(idxMinPeak:end));
        idxRefSample = find(shiftedWaveformAmps(idxMinPeak:end) ...
            >RATIO_VS_TALLEST.*ampMax, 1)+idxMinPeak-1;
        numSubFigs = 3;
    case 'FindPeaks'
        % TODO.
    otherwise
        error(['Peak searching method ',PEAK_SEARCH_METHOD,' undefined!']);
end

% Find the standard deviation and the mean for the specified part of the
% sample segment before the peak, and set the threshold accordingly.
idxRangeSampsToCons = floor(idxRefSample*relativeRangeForNoise(1)): ...
    ceil(idxRefSample*relativeRangeForNoise(2));
sigmaNoise = std(waveformToCons( idxRangeSampsToCons ));
meanNoise = mean(waveformToCons( idxRangeSampsToCons ));

if all(waveform>=0)
    % The amplitude case.
    thresholdMax = NUM_SIGMA_FOR_THRESHOLD*sigmaNoise + meanNoise;
    boolsEliminatedPts = abs(waveform)<thresholdMax;
else
    % The real / imaginary part case.
    thresholdMin = -NUM_SIGMA_FOR_THRESHOLD*sigmaNoise + meanNoise;
    thresholdMax = NUM_SIGMA_FOR_THRESHOLD*sigmaNoise + meanNoise;
    boolsEliminatedPts = (waveform<thresholdMax) & (waveform>thresholdMin);
end
% Set everything below the threshold to 0.
waveformThresholded = waveform;
waveformThresholded(boolsEliminatedPts) = 0;

%% Plots for debugging.
if flagDebug
    % Generate plots for debugging.
    hDebugFig = figure;
    subFigCounter = 1;
    
    if strcmp(PEAK_SEARCH_METHOD, 'TallEnoughAboveMedian')
        subplot(numSubFigs,1,subFigCounter); hold on;
        hShiftedWaveformAmp = plot((1:idxRefSample)./Fs, ...
            shiftedWaveformAmps(1:idxRefSample), '.', ...
            'Color', [1,1,1]*0.4);
        maxShifedWaveformAmp = plot([0 0;1 1].*idxRefSample./Fs, ...
            [1, -1;1 -1]*ampMax, 'r-.');
        hold off; grid minor; axis tight;
        legend([hShiftedWaveformAmp, maxShifedWaveformAmp(1)], ...
            'Used to compute delta', '+/- ampMax');
        transparentizeCurLegends;
        title('Samples for computing noise sigma');
        ylabel('Signal strength'); xlabel('time (s)');
        subFigCounter = subFigCounter+1;
    end
    
    subplot(numSubFigs,1,subFigCounter); hold on;
    hPtsNotCons = plot((1:idxRefSample)./Fs, ...
        waveform((1:idxRefSample)+numStartSampsToDiscard), '.', ...
        'Color', [1,1,1]*0.4);
    hPtsConsidered = plot(idxRangeSampsToCons./Fs, ...
        waveform(idxRangeSampsToCons+numStartSampsToDiscard), '*b');
    hold off; grid minor; axis tight;
    legend([hPtsConsidered, hPtsNotCons], ...
        'Used to compute delta', 'Other points before peak #1');
    transparentizeCurLegends;
    title('Samples for computing noise sigma');
    ylabel('Signal strength'); xlabel('time (s)');
    subFigCounter = subFigCounter+1;
    
    subplot(numSubFigs,1,subFigCounter); hold on;
    idxMaxToShow = min(floor(idxRefSample*2),length(waveform));
    waveformToShow = waveform(1:idxMaxToShow);
    plot((1:idxMaxToShow)./Fs, ...
        waveformToShow, '.', ...
        'Color', 'b');
    indicesDiscared = find(abs(waveformToShow)<thresholdMax);
    hEliminatedPts = plot(indicesDiscared./Fs, ...
        waveformToShow(indicesDiscared), '.', ...
        'Color', [1,1,1]*0.5);
    hTwiceMedian = ...
        plot([0 0;1 1].*(2*idxRefSample+1)./Fs, ...
        [1, -1;1 -1]*twiceMedianShifted, 'r-.');
    grid minor; axis tight; finalAxis = axis;
    hDiscardRegion = ...
        patch([0,0,2*idxRefSample+1,2*idxRefSample+1]./Fs, ...
        [-thresholdMax,thresholdMax,thresholdMax,-thresholdMax], ...
        [1,1,1]*0.9);
    set(hDiscardRegion, 'LineStyle', 'None');
    uistack(hDiscardRegion,'bottom')
    hold off;
    axis(finalAxis);
    legend([hDiscardRegion, hEliminatedPts, hTwiceMedian(1)], ...
        ['Elimination Region (', num2str(NUM_SIGMA_FOR_THRESHOLD),'*sigma+mean)'], 'Eliminated Samples', '+/- twiceMedianShifted');
    transparentizeCurLegends;
    title('Eliminated samples (Set to 0)');
    ylabel('Signal strength'); xlabel('time (s)');
    subFigCounter = subFigCounter+1;
end

end
% EOF