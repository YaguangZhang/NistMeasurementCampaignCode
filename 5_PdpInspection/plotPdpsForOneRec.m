function [ hPdpFig, nsToPlot, signalAmp, lowPassedSig, indexRangeShown, ...
    hNoiseEliDebugFig] ...
    = plotPdpsForOneRec(sigOutFile, F_S, segmentRange, ...
    flagGenerateNoiseEliDebugFig)
%PLOTPDPSFORONEREC Plot the PDP overview plot for one signal recording
%file.
%
% Note: for applying noise elimination before showing the signal,
% parameters like Fs and NUM_SIGMA_FOR_THRESHOLD need to be defined in the
% base workspace.
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
%       - flagGenerateNoiseEliDebugFig
%         Set this to be true if it is necessary to generate the debug plot
%         for noise elimination.
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
%       - hPdpFig
%         The handler to the figure generated.
%       - nsToPlot, signalAmp
%         Time points in ns and the cooresponding signal sample amplitudes
%         in the plot. If everything works as expected, these should
%         capture the tallest peak in the signal.
%       - lowPassedSig
%         The original signal, cooresponding to signalAmp, before noise
%         elimination and amplitude computation.
%       - indexRangeShown
%         The sample index range for PDP segment captured in the original
%         signal file.
%       - hNoiseEliDebugFig
%         The handler for the noise elimination debug figure.
%
% Update 06/28/2018: Added a new output (indexRangeShown) for debugging.
%
% Update 06/22/2018: Added a new output (lowPassedSig) for debugging.
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

if ~exist('flagGenerateNoiseEliDebugFig', 'var')
    flagGenerateNoiseEliDebugFig = false;
end

if (exist('segmentRange', 'var') && (~any(isnan(segmentRange))))
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
% Output for extra visualization if necessary.
lowPassedSig = curSignal;

% Noise elimination.

% noiseEliminationFct = @(waveform) thresholdWaveform(abs(waveform), ...
%     flagGenerateNoiseEliDebugFig);
% [~, boolsEliminatedPts, hNoiseEliDebugFig] =
% noiseEliminationFct(curSignal);
[~, boolsEliminatedPts, hNoiseEliDebugFig] ...
    = thresholdWaveform(abs(curSignal), ...
    flagGenerateNoiseEliDebugFig);
curSignalEliminated = curSignal;
curSignalEliminated(boolsEliminatedPts) = 0;

% Also get rid of everything below the USRP noise floor if
% USRP_NOISE_FLOOR_V is specified in the base workspace.
if evalin('base','exist(''USRP_NOISE_FLOOR_V'', ''var'')')    
    USRP_NOISE_FLOOR_V = evalin('base', 'USRP_NOISE_FLOOR_V');
    %     disp(['USRP_NOISE_FLOOR_V is set to be ', ...
    %         num2str(USRP_NOISE_FLOOR_V)])
    curSignalEliminated...
        (abs(curSignalEliminated)<USRP_NOISE_FLOOR_V) = 0;
else
    disp('USRP_NOISE_FLOOR_V is not defined.')
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
        signalToShow = signalToShow(end:-1:1);
        lowPassedSig = lowPassedSig(end:-1:1);
    end
else
    disp('FLAG_PDP_TIME_REVERSED is not defined.')
end

if evalin('base','exist(''SLIDE_FACTOR'', ''var'')')
    SLIDE_FACTOR = evalin('base', 'SLIDE_FACTOR');
    [hPdpFig, nsToPlot, signalAmp, signalIdxRange] ...
        = plotOnePresentSignalAmp(signalToShow, ...
        numPreSamples, numPostSamples, figureSupTitle, F_S, SLIDE_FACTOR);
else
    disp('SLIDE_FACTOR is not defined.')
    [hPdpFig, nsToPlot, signalAmp, signalIdxRange] ...
        = plotOnePresentSignalAmp(signalToShow, ...
        numPreSamples, numPostSamples, figureSupTitle, F_S);
end
if isempty(signalAmp)
    lowPassedSig = signalAmp;
else
    lowPassedSig = lowPassedSig(signalIdxRange(1):signalIdxRange(2));
end

rangeIndices = range(1):range(2);
if ~isempty(signalIdxRange)
    indexRangeShown = rangeIndices([signalIdxRange(1), signalIdxRange(2)]);
else
    indexRangeShown = [];
end
end
% EOF