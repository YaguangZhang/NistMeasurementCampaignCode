function [ hPsdInDb, hPsdInDbZoomedIn ] ...
    = plotPsdInDbWithLPF( powerSpectralDen, f, maxFreqPassed, ...
    minFreqPassed, estimatedSnr )
%PLOTPSDINDBWITHLPF Plot the PSD (in dB) for the input power spectral
%density with frequency range f.
%
% And maxFreqPassed is the cutoff frequency of the input LPF.
% Correspondingly, the signal frequency range that can pass will also be
% plotted.
%
% Update 20170919: Also add an HPF to remove the very low-frequency / DC
% components.
%
% Optionally, if estimatedSnr is specified, it will be shown on the plot,
% too.
%
% Yaguang Zhang, Purdue, 09/13/2017

hPsdInDb = figure; hold on;
powerSpectralDenIndB = 10*log10(powerSpectralDen);
hPowerSpectralDenIndB = plot(f, powerSpectralDenIndB);
curAxis = axis; curAxis(1) = f(1); curAxis(2) = f(end);
if ~isinf(maxFreqPassed)
    hLPFIndB = plot([maxFreqPassed, -maxFreqPassed; ...
        maxFreqPassed, -maxFreqPassed], ...
        [curAxis(3),curAxis(3); ...
        curAxis(4),curAxis(4)], '-.r');
    x = [maxFreqPassed maxFreqPassed f(end) f(end)];
    y = [curAxis(3) curAxis(4) curAxis(4) curAxis(3)];
    patch(x,y,[1,1,1].*0.6,'FaceAlpha',0.3,'LineStyle','none');
    x = [-maxFreqPassed -maxFreqPassed f(1) f(1)];
    y = [curAxis(3) curAxis(4) curAxis(4) curAxis(3)];
    patch(x,y,[1,1,1].*0.6,'FaceAlpha',0.3,'LineStyle','none');
end
if minFreqPassed>0
    hHPFIndB = plot([minFreqPassed, -minFreqPassed; ...
        minFreqPassed, -minFreqPassed], ...
        [curAxis(3),curAxis(3); ...
        curAxis(4),curAxis(4)], '-.k');
    x = [-minFreqPassed -minFreqPassed ...
        minFreqPassed minFreqPassed];
    y = [curAxis(3) curAxis(4) curAxis(4) curAxis(3)];
    patch(x,y,[1,1,1].*0.6,'FaceAlpha',0.3,'LineStyle','none');
end
if nargin>4
    text(min(maxFreqPassed, 50000), mean([curAxis(3) curAxis(4)]), ...
        ['Estimated SNR = ', ...
        num2str(estimatedSnr, '%.2f')]);
end
hold off;
xlabel('f (Hz)'); ylabel('Estimated PSD (V^2/Hz in dB)');
axis(curAxis);
if ~isinf(maxFreqPassed) && minFreqPassed>0
    curLegends = legend([hPowerSpectralDenIndB, hLPFIndB(1), hHPFIndB(1)], ...
        'PSD (dB)', 'LPF', 'HPF');
elseif isinf(maxFreqPassed) && minFreqPassed>0
    curLegends = legend([hPowerSpectralDenIndB, hHPFIndB(1)], ...
        'PSD (dB)', 'HPF');
elseif ~isinf(maxFreqPassed) && minFreqPassed<=0
    curLegends = legend([hPowerSpectralDenIndB, hLPFIndB(1)], ...
        'PSD (dB)', 'LPF');
else
    legend(hPowerSpectralDenIndB, 'PSD (dB)');
end
transparentizeCurLegends;
grid minor;

% Generate a zoomed-in version of the plot when necessary.
if nargout>1
    hPsdInDbZoomedIn = cloneFig(hPsdInDb);
    % Need to fix the legends.
    set(findobj(hPsdInDbZoomedIn, 'Type', 'legend'), 'String', ...
        get(curLegends, 'String'));
    curAxis = axis;
    axis([[-minFreqPassed minFreqPassed].*6 curAxis(3:4)]);
end

end
% EOF