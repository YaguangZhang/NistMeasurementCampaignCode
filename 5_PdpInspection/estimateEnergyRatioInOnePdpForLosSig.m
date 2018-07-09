function [ energyRatioForLosSig ] ...
    = estimateEnergyRatioInOnePdpForLosSig(...
    timesForOnePdp, samAmpsForOnePdp, ...
    fullPathToSavePlot, lowPassedSigForOnePdp)
%ESTIMATEENERGYRATIOFORLOSSIG Estimate the energy ratio of the LOS signal
%(the first arrival peak) for one PDP.
%
%   Inputs:
%       - timesForOnePdp, samAmpsForOnePdp
%         The time points and the corresponding sample amplitudes for the
%         input PDP. Note that we assume that the PDP samples have been
%         gone through LPF and noise eliminiation already.
%       - fullPathToSavePlot
%         Optional. The full path, including the file name, to save a debug
%         figure. If absent, no figure will be generated.
%       - lowPassedSigForOnePdp
%
%
%   Output:
%       - energyRatioForLosSig
%         A scalar in [0,1] to present the estimated signal energy ratio of
%         the first peak.
%
%   We will find the peaks in the PDP and then estimate each peak's energy
%   simply by (peak amplitude)^2.
%
%   Update 20180622: also show the oringinal singal amplitude for
%   debugging.
%
% Yaguang Zhang, Purdue, 06/06/2018

% If one peak is too close to the previous one, and they are of similar
% height, we will ignore it.
MIN_NUM_SAMPS_BETWEEN_VALID_PEAKS = 10;
MAX_HEIGHT_PERCENT_DIFF_OK_TO_IGNORE = 0.01;

numSamAmpsForOnePdp = length(samAmpsForOnePdp);

if isempty(samAmpsForOnePdp)
    energyRatioForLosSig=nan;
    [pks, locs] = deal([]);
    boolsValidSigPeaks = true(1, length(pks));
else
    [pks, locs] = findpeaks(samAmpsForOnePdp);
    boolsValidSigPeaks = true(1, length(pks));
    % The LoS peak should be at least 20% of the highest signal received.
    idxLoSPeak = find(pks./max(pks)>=0.2, 1);
    boolsValidSigPeaks(1:(idxLoSPeak-1)) = false;
    
    locPreValidPeak = locs(idxLoSPeak);
    for idxPeak = (idxLoSPeak+1):length(pks)
        if boolsValidSigPeaks(idxPeak)
            curLoc = locs(idxPeak);
            samIdxPrev = curLoc-1;
            samIdxPost = curLoc+1;
            % Peaks after LoS should not be isolated. We get rid of samples
            % next to any zero samples.
            if ((samIdxPrev>0 && samAmpsForOnePdp(samIdxPrev)==0) ...
                    || (samIdxPost<=numSamAmpsForOnePdp ...
                    && samAmpsForOnePdp(samIdxPost)==0))
                boolsValidSigPeaks(idxPeak) = false;
            end
            % We also get rid of peaks too close to each other.
            if ((curLoc-locPreValidPeak)...
                    <MIN_NUM_SAMPS_BETWEEN_VALID_PEAKS) ...
                    && (abs(samAmpsForOnePdp(curLoc) ...
                    -samAmpsForOnePdp(locPreValidPeak)) ...
                    ./min(abs(samAmpsForOnePdp(curLoc)), ...
                    abs(samAmpsForOnePdp(locPreValidPeak)))...
                    <MAX_HEIGHT_PERCENT_DIFF_OK_TO_IGNORE)
                    boolsValidSigPeaks(idxPeak) = false;
            else
                locPreValidPeak = curLoc;
            end
        end
    end
    
    % The LoS peak should be the first meaningful signal received.
    energyRatioForLosSig = (pks(idxLoSPeak)).^2 ...
        /sum(pks(boolsValidSigPeaks).^2);
end

if exist('fullPathToSavePlot', 'var')
    % Generate a plot for the zoomed-in version of the PDP and the
    % estimated energy ratio for the LOS signal labeld.
    extraNumSampsPerSide = 1000; % Only plot a little more samples.
    indicesSampsToShow = max([1, min(locs)-extraNumSampsPerSide]): ...
        min([max(locs)+extraNumSampsPerSide, numSamAmpsForOnePdp]);
    timesToShow = timesForOnePdp(indicesSampsToShow);
    samAmpsToShow = samAmpsForOnePdp(indicesSampsToShow);
    
    timeLocs = timesForOnePdp(locs);
    
    % Plot.
    hZoomedInPdp = figure; hold on;
    % Only plot valid signal peaks.
    hPeaks = plot(timeLocs(boolsValidSigPeaks), ...
        pks(boolsValidSigPeaks), 'ro', 'LineWidth', 1);
    
    % Double the x range if there are multiple peaks.
    %     if sum(boolsValidSigPeaks)>1
    %         curAxis = axis; curAxis(1:2) = curAxis(1:2) ...
    %             + [-1, 1].*0.5.*(curAxis(2)-curAxis(1));
    %     end
    
    hAmp = plot(timesToShow, samAmpsToShow, 'b.');
    % Set x range according to extraNumSampsPerSide if there is not enough
    % peaks.
    %     if sum(boolsValidSigPeaks)<=1
    curAxis = axis;
    %     end
    
    if exist('lowPassedSigForOnePdp', 'var')
        hOri = plot(timesForOnePdp, abs(lowPassedSigForOnePdp), '.', ...
            'Color', ones(1,3)*0.7);
        legend([hAmp, hPeaks, hOri], ...
            'Sample amplitude', 'Peaks', 'Eliminated');
    else
        legend([hAmp, hPeaks], 'Sample amplitude', 'Peaks');
    end
    
    title({'Peaks Detected (Amplitude in Volt)'; ...
        ['Energy ratio for the first peak = ', ...
        num2str(energyRatioForLosSig, '%.4f')]});
    hold off; axis tight; transparentizeCurLegends;
    grid on; xlabel('Time (ms)'); ylabel('Amplitude (Volage)');
    % Focus on the peaks.
    axis([curAxis(1:2) 0 curAxis(4)]);
    uistack(hAmp, 'top'); uistack(hPeaks, 'top')
    
    saveas(hZoomedInPdp, fullPathToSavePlot);
    close(hZoomedInPdp);
end

end
% EOF