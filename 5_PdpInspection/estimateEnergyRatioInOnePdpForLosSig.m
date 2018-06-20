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
    idxLoSPeak = 1;
else
    [pks, locs] = findpeaks(samAmpsForOnePdp);
    % The LoS peak should be at least 20% of the highest signal received.
    idxLoSPeak = find(pks./max(pks)>=0.2, 1);
    % The LoS peak should be the first meaningful signal received.    
    energyRatioForLosSig = (pks(idxLoSPeak)).^2 ...
        /sum(pks(idxLoSPeak:end).^2);
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
    hAmp = plot(timesToShow, samAmpsToShow, 'b.');
    % Hide peaks before the LoS one.
    hPeaks = plot(timeLocs(idxLoSPeak:end), pks(idxLoSPeak:end), 'ro');
    legend([hAmp, hPeaks], 'Sample amplitude', 'Peaks');
    title({'Peaks Detected (Amplitude in Volt)'; ...
        ['Energy ratio for the first peak = ', ...
        num2str(energyRatioForLosSig, '%.4f')]});
    hold off; axis tight; transparentizeCurLegends;
    grid on; xlabel('Time (ms)'); ylabel('Amplitude (Volage)');
    
    saveas(hZoomedInPdp, fullPathToSavePlot);
    close(hZoomedInPdp);
end

end
% EOF