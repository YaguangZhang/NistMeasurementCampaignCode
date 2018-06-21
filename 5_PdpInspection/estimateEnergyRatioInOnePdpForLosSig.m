function [ energyRatioForLosSig ] ...
    = estimateEnergyRatioInOnePdpForLosSig(...
    timesForOnePdp, samAmpsForOnePdp, ...
    fullPathToSavePlot)
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
%
%   Output:
%       - energyRatioForLosSig
%         A scalar in [0,1] to present the estimated signal energy ratio of
%         the first peak.
%
%   We will find the peaks in the PDP and then estimate each peak's energy
%   simply by (peak amplitude)^2.
%
% Yaguang Zhang, Purdue, 06/06/2018

if isempty(samAmpsForOnePdp)
    energyRatioForLosSig=nan;
    [pks, locs] = deal([]);
    boolsValidSigPeaks = ones(1, length(pks));
else
    [pks, locs] = findpeaks(samAmpsForOnePdp);
    % The LoS peak should be at least 20% of the highest signal received.
    boolsValidSigPeaks = pks./max(pks)>=0.2;
    idxLoSPeak = find(boolsValidSigPeaks, 1);
    % The LoS peak should be the first meaningful signal received.    
    energyRatioForLosSig = (pks(idxLoSPeak)).^2 ...
        /sum(pks(boolsValidSigPeaks).^2);
end

if exist('fullPathToSavePlot', 'var')
    % Generate a plot for the zoomed-in version of the PDP and the
    % estimated energy ratio for the LOS signal labeld.
    extraNumSampsPerSide = 10;
    indicesSampsToShow = max([1, min(locs)-extraNumSampsPerSide]): ...
        min([max(locs)+extraNumSampsPerSide, length(samAmpsForOnePdp)]);
    timesToShow = timesForOnePdp(indicesSampsToShow);
    samAmpsToShow = samAmpsForOnePdp(indicesSampsToShow);
    
    timeLocs = timesForOnePdp(locs);
    
    % Plot.
    hZoomedInPdp = figure; hold on;
    % Only plot valid signal peaks.
    hPeaks = plot(timeLocs(boolsValidSigPeaks), ...
        pks(boolsValidSigPeaks), 'ro', 'LineWidth', 1);
    curAxis = axis;
    hAmp = plot(timesToShow, samAmpsToShow, 'b.');    
    legend([hAmp, hPeaks], 'Sample amplitude', 'Peaks');
    title({'Peaks Detected (Amplitude in Volt)'; ...
        ['Energy ratio for the first peak = ', ...
        num2str(energyRatioForLosSig, '%.4f')]});
    hold off; axis tight; transparentizeCurLegends;
    grid on; xlabel('Time (ms)'); ylabel('Amplitude (Volage)');
    % Focus on the peaks.
    axis([curAxis(1:2) 0 curAxis(4)]);
    uistack(hPeaks, 'top')
    
    saveas(hZoomedInPdp, fullPathToSavePlot);
    close(hZoomedInPdp);
end

end
% EOF