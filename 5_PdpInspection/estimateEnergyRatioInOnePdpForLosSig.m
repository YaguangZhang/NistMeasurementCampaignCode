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

boolsValidSigPeaks = true(1, length(pks));
if isempty(samAmpsForOnePdp)
    energyRatioForLosSig=nan;
    [pks, locs] = deal([]);
else
    [pks, locs] = findpeaks(samAmpsForOnePdp);
    % The LoS peak should be at least 20% of the highest signal received.
    idxLoSPeak = find(pks./max(pks)>=0.2, 1);
    boolsValidSigPeaks(1:(idxLoSPeak-1)) = false;
    % Peaks after LoS should not be isolated. We get rid of samples next to
    % any zero samples.
    numSamAmpsForOnePdp = length(samAmpsForOnePdp);
    for idxPeak = (idxLoSPeak+1):length(pks)
        if boolsValidSigPeaks(idxPeak)
            samIdxPrev = locs(idxPeak)-1;
            samIdxPost = locs(idxPeak)+1;
            if ((samIdxPrev>0 && samAmpsForOnePdp(samIdxPrev)==0) ...
                    || (samIdxPost<=numSamAmpsForOnePdp ...
                    && samAmpsForOnePdp(samIdxPost)==0))
                boolsValidSigPeaks(idxPeak) = false;
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
    extraNumSampsPerSide = 500; % Only plot a little more samples.
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
    if sum(boolsValidSigPeaks)>1
        curAxis = axis;
        curAxis(1:2) = curAxis(1:2) ...
            + [-1, 1].*0.5.*(curAxis(2)-curAxis(1));
    end
    
    hAmp = plot(timesToShow, samAmpsToShow, 'b.');
    % Set x range according to extraNumSampsPerSide if there is not enough
    % peaks.
    if sum(boolsValidSigPeaks)<=1
        curAxis = axis;
    end
    
    hOri = plot(timesToShow, abs(lowPassedSigForOnePdp), '.', ...
        'Color', ones(1,3)*0.7);
    
    legend([hAmp, hPeaks, hOri], 'Sample amplitude', 'Peaks', 'Eliminated');
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