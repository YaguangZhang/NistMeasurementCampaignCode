function [ hFig, nsToPlot, signalAmp, signalIdxRange ] ...
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
%     Optional. A string to specify the name of the figure.
%   - Fs
%     Optional. A number to specify the sample rate for the input signal.
%     Used to properly set the x axis / time line labels. If not present, a
%     default value of 1.04 MSamples/s will be used.
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
% Update 20180622: Also output the index range for the signal to be
% plotted.
%
% Yaguang Zhang, Purdue, 04/30/2018

% By default, find the flag flagGenFigSilently in the base workspace. If
% not found, flagGenFigSilently will be set to false.
try
    flagGenFigSilently = evalin('base', 'flagGenFigSilently');
catch
    flagGenFigSilently = false;
end

FLAG_SUBTITLES = true;

if nargin<5
    Fs = 1.04*10^6;
end

RATIO_VS_TALLEST = 0.9;

signalAmp = abs(signal);

[valueMax, ~] = max(signalAmp);
% Find the first peak that is tall enough.
idxRefSample = find(signalAmp>RATIO_VS_TALLEST.*valueMax, 1);

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

if exist('figureName', 'var') && isstring(figureName)
    hFig = figure('Name',figureName, 'visible', ~flagGenFigSilently);
else
    hFig = figure('visible', ~flagGenFigSilently);
end

% Convert sample number to ns.
nsToPlot = (1:(maxIdxToPlot-minIdxToPlot+1))./Fs.*(10.^9);
if nargin>5
    nsToPlot = nsToPlot./slideFactor;
end

% Output the index range of the signal of interest.
signalIdxRange = [minIdxToPlot, maxIdxToPlot];

subplot(1,1,1); hold on;
signalAmp = abs(signal(minIdxToPlot:maxIdxToPlot));
hAmp = plot(nsToPlot, signalAmp, 'b.');
if FLAG_SUBTITLES
    title('First Detected Signal (Amplitude in Volt)')
end
hold off; legend(hAmp, 'Amplitude'); axis tight; grid on;
if nargin>3
    suptitle(figureName);
end
xlabel('Time (ns)');

set(0,'DefaultTextInterpreter',curDefulatTextInt);
% EOF