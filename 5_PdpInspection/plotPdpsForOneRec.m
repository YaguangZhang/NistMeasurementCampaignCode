function [ hFig ] = plotPdpsForOneRec(sigOutFile, F_S, segmentRange)
%PLOTPDPSFORONEREC Plot the PDP overview plot for one signal recording
%file.
%
%   Inputs:
%       - sigOutFile
%         One struct of the output array from rdir specifying where the
%         signal .out file is.
%       - F_S
%         The GnuRadio sample rate for the singal recordings. The x axis
%         will be converted from sample # to time accordingly.
%       - segmentRange
%         Optional. [idxStart, idxEnd] for what segment to look at
%         (including the start and end samples). If not specified, the
%         whole recording will be looked at.
%
% Update 05/07/2018: Added support for segment of the signal.
%
% Yaguang Zhang, Purdue, 04/30/2018

if exist('segmentRange', 'var')
    % We will load in ~1 second, if possible, of the signal recording
    % segment from the center.
    numSam = segmentRange(2) - segmentRange(1) + 1;
    centerSam = segmentRange(1)+floor(numSam/2);
else
    % We will load in ~1 second, if possible, of the whoe signal recording
    % from the center.
    numSam = countComplexBinary(sigOutFile.name); % Total number of samples.
    centerSam = floor(numSam/2);
end
% Number of samples to load.
countSam = min(F_S, numSam); % At most 1s of the signal.

shiftSam = floor(countSam/2);
range = [centerSam-shiftSam, centerSam+shiftSam];

curSignal = readComplexBinaryInRange (sigOutFile.name, range);
[~, figureSupTitle, ~] = fileparts(sigOutFile.name);

% For constructing the LPF, which will be applied to the signal loaded
% before any other operations.
Fp  = 10e3;   	% 10 kHz passband-edge frequency
Fst = 12e3;     % Transition Width = Fst - FpmaxFreqPassed
Ap = 0.01;      % Allowed peak-to-peak ripple
Ast = 80;       % Stopband attenuation
% Filter the signal with a LPF.
lpfComplex = dsp.LowpassFilter('SampleRate', F_S, ...
    'FilterType', 'FIR', 'PassbandFrequency', Fp, ...
    'StopbandFrequency', Fst, ...
    'PassbandRipple', Ap, ...
    'StopbandAttenuation', Ast ...
    );
release(lpfComplex);
curSignal = lpfComplex(curSignal);

% Plot the signals. We will try to find the "tallest" bump for each
% measurement.
numPreSamples = floor(F_S*0.02); % 0.02s
numPostSamples = floor(F_S*0.04); % 0.04s


if evalin('base','exist(''SLIDE_FACTOR'', ''var'')')
    SLIDE_FACTOR = evalin('base', 'SLIDE_FACTOR');
    hFig = plotOnePresentSignalAmp(curSignal, ...
        numPreSamples, numPostSamples, figureSupTitle, F_S, SLIDE_FACTOR);
else
    hFig = plotOnePresentSignalAmp(curSignal, ...
        numPreSamples, numPostSamples, figureSupTitle, F_S);
end
transparentizeCurLegends; grid on; xlabel('Time (ms)');

end
% EOF