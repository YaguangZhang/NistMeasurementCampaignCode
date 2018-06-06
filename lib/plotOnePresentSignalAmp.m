function [ hFig, msToPlot, signalAmp ] ...
    = plotOnePresentSignalAmp( signal, ...
    numPreSamples, numPostSamples, figureName, Fs, slideFactor)
%PLOTONEPRESENTSIGNALAMP Plot the tallest bump of the signal in amplitude.
%
% Inputs:
%   - signal
%     A vector of complex numbers.
%   - numPreSamples, numPostSamples
%     Integers. The range of signal to be plotted can be adjusted using
%     numPreSamples and numPostSamples. The samples before and after the
%     first tallest signal sample found will be plotted accordingly if
%     there are enough samples available.
%   - figureName
%     Optional. A string to specify the figure's name.
%   - Fs
%     Optional. A number to specify the sample rate for the input signal.
%     Used to properly set the x axis / time line labels. Default
%     1.04MSamples/s.
%
% Update 20180518: Considered the time dilation effect when setting the x
% labels.
%
% New inputs:
%   - slideFactor
%     (Optional) The slider factor for the sliding correlator channel
%     souder used. This is determined by the psuedo-noise sequence
%     generator clock frequencies for the TX and the RX. When this is
%     present, the time line (x axis) will be adjusted according to
%     compensate for the slider factor.
%
% Yaguang Zhang, Purdue, 04/30/2018

FLAG_SUBTITLES = true;

if nargin<5
    Fs = 1.04*10^6;
end

RATIO_VS_TALLEST = 0.9;

signalReal = real(signal);
signalImag = imag(signal);

[valueMax, ~] = max(real(signal));
% Find the first peak that is tall enough.
idxRefSample = find(signalReal>RATIO_VS_TALLEST.*valueMax, 1);

minIdxToPlot = idxRefSample - numPreSamples;
if(minIdxToPlot<1)
    minIdxToPlot = 1;
end
maxIdxToPlot = idxRefSample + numPostSamples;
if(maxIdxToPlot>length(signal))
    maxIdxToPlot = length(signal);
end

% Disable the interpreter temperorily.
curDefulatTextInt = get(0,'DefaultTextInterpreter');
set(0,'DefaultTextInterpreter','none');

if nargin>3
    hFig = figure('Name',figureName);
else
    hFig = figure;
end

% Convert sample number to ms.
msToPlot = (1:(maxIdxToPlot-minIdxToPlot+1))./Fs.*1000;
if nargin>5
    msToPlot = msToPlot./slideFactor;
end

subplot(1,1,1); hold on;
signalAmp = abs(signal(minIdxToPlot:maxIdxToPlot));
hAmp = plot(msToPlot, signalAmp, 'b.');
if FLAG_SUBTITLES
    title('First Detected Signal (Amplitude in Volt)')
end
hold off; legend(hAmp, 'Amplitude'); axis tight;

if nargin>3
    suptitle(figureName);
end

set(0,'DefaultTextInterpreter',curDefulatTextInt);
% EOF