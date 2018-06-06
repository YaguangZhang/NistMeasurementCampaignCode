function [ hFig, msToPlot, signalAmp ] ...
    = plotPdpsForOneRec(sigOutFile, F_S, segmentRange)
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
%   Optional parameters in the base workspace:
%       - FLAG_PDP_TIME_REVERSED
%         If present and set to be true, the PDP plot will be generated
%         with time reversed back.
%       - SLIDE_FACTOR
%         If present, the value will be used to adjust the labels for the
%         time (x) axis.
%
%   Outputs:
%       - hFig
%         The handler to the figure generated.
%       - msToPlot, signalAmp
%         Time points in ms and the cooresponding signal sample amplitudes
%         in the plot. If everything works as expected, these should
%         capture the tallest peak in the signal.
%
% Update 06/06/2018: Added noise eliminate before plotting.
%
% Update 05/31/2018: Added extra plotting functions controled by optional
% parameters in the base workspace.
%
% Update 05/07/2018: Added support for segment of the signal.
%
% Yaguang Zhang, Purdue, 04/30/2018

% Add path to thresholdWaveForm.m.
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '1_Calibration'));

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

% Noise elimination.
noiseEliminationFct = @(waveform) thresholdWaveform(abs(waveform));
[~, boolsEliminatedPts] = noiseEliminationFct(curSignal);
curSignalEliminated = curSignal;
curSignalEliminated(boolsEliminatedPts) = 0;

% Also get rid of everything below the USRP noise floor if
% USRP_NOISE_FLOOR_V is specified in the base workspace.
if evalin('base','exist(''USRP_NOISE_FLOOR_V'', ''var'')')
    USRP_NOISE_FLOOR_V = evalin('base', 'USRP_NOISE_FLOOR_V');
    curSignalEliminated...
        (abs(curSignalEliminated)<USRP_NOISE_FLOOR_V) = 0;
end

% Plot the signals. We will try to find the "tallest" bump for each
% measurement.
numPreSamples = floor(F_S*0.02); % 0.02s
numPostSamples = floor(F_S*0.04); % 0.04s

% Reverse the signal if necessary.
signalToShow = curSignalEliminated;
if evalin('base','exist(''FLAG_PDP_TIME_REVERSED'', ''var'')')
    FLAG_PDP_TIME_REVERSED = evalin('base', 'FLAG_PDP_TIME_REVERSED');
    if FLAG_PDP_TIME_REVERSED
        signalToShow = curSignalEliminated(end:-1:1);
    end
end

if evalin('base','exist(''SLIDE_FACTOR'', ''var'')')
    SLIDE_FACTOR = evalin('base', 'SLIDE_FACTOR');
    [hFig, msToPlot, signalAmp] = plotOnePresentSignalAmp(signalToShow, ...
        numPreSamples, numPostSamples, figureSupTitle, F_S, SLIDE_FACTOR);
else
    [hFig, msToPlot, signalAmp] = plotOnePresentSignalAmp(signalToShow, ...
        numPreSamples, numPostSamples, figureSupTitle, F_S);
end
transparentizeCurLegends; grid on; xlabel('Time (ms)');

end
% EOF