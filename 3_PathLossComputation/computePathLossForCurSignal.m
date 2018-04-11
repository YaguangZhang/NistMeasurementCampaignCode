function [ pathLossInDb ] ...
    = computePathLossForCurSignal(curSignal, txPower, ...
    rxGain, noiseEliminationFct, powerShiftsForCali, FlagCutHead)
%COMPUTEPATHLOSSFORCURSIGNAL Compute the path loss for the input complex
%array curSignal. Note that here antGain = 0, i.e. the path actually
%includes the TX & RX antennas, because actual TX and RX ant. gains will
%depend on the TX & RX setups.
%
% Inputs:
%   - curSignal
%     A complex array representing the Rx signal.
%   - txPower
%     TX power (after upconverter) in dBm.
%   - rxGain
%     A scalar. The Gnu Radio gain (in dB) for curSignal.
%   - noiseEliminationFct
%     A function to run the loaded signal through for a noise-eliminated
%     copy of the signal. The syntax of it should be like:
%         [signalNoiseEliminated, boolsEliminatedPts] = ...
%             noiseEliminationFct(signalComplexArray)
%   - powerShiftsForCali
%     Essentially the parameter b in y=ax+b for the calibration line
%     corresponding to the current rxGain.
%   - FlagCutHead
%     Optional. True by default. Set this to false if there is no need to
%     discard the heading numStartSampsToDiscard samples.
%
% Procedures below will be carried out one by one:
%    (1) LPF
%        A pre-filtering procedure to reduce noise in frequency domain.
%    (2) Shrink the input signal sequence
%        Only use a segment of the signal for path loss computation.
%    (3) Noise elimination
%        Further reduce noise in time domain.
%    (4) Power calculation
%        Calculate the power vis PSD.
%    (5) Rx calibration (6) Antena normalization
%
% Update 04/10/2018: Adjusted parameters for the NIST dataset.
%
% Yaguang Zhang, Purdue, 04/09/2018

%% Parameters

% By default, we need to discard the first numStartSampsToDiscard samples
% to avoid the warm-up stage of USRP. Note that this is not necessary for
% most of the Conti measurements.
if nargin < 6
    FlagCutHead = true;
end

% For constructing the LPF, which will be applied to the signal loaded
% before any other operations.
Fp  = 60e3;   	% 60 kHz passband-edge frequency
Fst = 65e3;     % Transition Width = Fst - FpmaxFreqPassed
Ap = 0.01;      % Allowed peak-to-peak ripple
Ast = 80;       % Stopband attenuation

% Low pass filter for the PSD. Tried before: 46000; 39500.
maxFreqPassed = 10000; % In Hz; 10kHz for NIST data.
% High pass filter to remove the DC component.
minFreqPassed = 1; % In Hz.

% After discarding these samples, furthermore only keep the middle part of
% the signal for calibration.
timeLengthAtCenterToUse = 1; % In second.

% Sample rate used for GnuRadio.
try
    Fs = evalin('base', 'Fs');
catch
    warning('GnuRadio sample frequency Fs not found in the base workspace.')
    warning('Will use the default value 1.04 * 10^6.')
    Fs = 2 * 10^6;
end

% Number of samples to discard at the beginning.
numStartSampsToDiscard = 0.1*Fs; % ~0.1s

% The downconverter gain at the RX side.
try
    downConvGainInDb = evalin('base', 'DOWNCONVERTER_GAIN_IN_DB');
catch
    warning('The gain for the downconverter DOWNCONVERTER_GAIN_IN_DB not found in the base workspace.')
    warning('Will set it to 0.')
    downConvGainInDb = 0;
end

%% LPF

% Pre-filter the input with a LPF.
lpfComplex = dsp.LowpassFilter('SampleRate', Fs, ...
    'FilterType', 'FIR', 'PassbandFrequency', Fp, ...
    'StopbandFrequency', Fst, ...
    'PassbandRipple', Ap, ...
    'StopbandAttenuation', Ast ...
    );
release(lpfComplex);
curSignal = lpfComplex(curSignal);

%% Get a Segment of the Signal for Path Loss Computation

if FlagCutHead
    % Discard the first numStartSampsToDiscard of samples.
    curSignal = curSignal((numStartSampsToDiscard+1):end);
end
% Further more, only keep the middle part for calibration.
numSampsToKeep = ceil(timeLengthAtCenterToUse*Fs);
numSampsCurSeries = length(curSignal);
if numSampsToKeep > numSampsCurSeries
    warning('There are not enough samples to keep. We will use all remaining ones.');
else
    idxRangeToKeep = floor(0.5.*numSampsCurSeries ...
        + [-1,1].*numSampsToKeep./2);
    curSignal = curSignal(max(idxRangeToKeep(1),1) ...
        :min(idxRangeToKeep(2), numSampsCurSeries));
end
% Make sure we end up with even number of samples.
if mod(length(curSignal),2)==1
    curSignal = curSignal(1:(end-1));
end

%% Noise Elimination

% Noise elimination.
[~, boolsEliminatedPts] = ...
    noiseEliminationFct(curSignal);
curSignalEliminated = curSignal;
curSignalEliminated(boolsEliminatedPts) = 0;

%% Calculate Power

% For the signal to process, i.e. the noise eliminiated signal, compute the
% PSD.
X = curSignalEliminated;
L = length(X);
% FFT results.
Y = fftshift(fft(X));
% Frequency domain.
f = (-L/2:L/2-1)*(Fs/L);
idxDC = L/2+1;
% PSD.
powerSpectralDen = abs(Y).^2/L;

% Compute the power.
boolsFPassed = abs(f)<=maxFreqPassed ...
    & abs(f)>=minFreqPassed;
% Compute the power by integral. Note that we will always discard the DC
% component here (although it may be passed by the filters).
psdPassed = powerSpectralDen;
psdPassed(~boolsFPassed) = 0;
psdPassed(idxDC) = 0;
calcP = trapz(f, psdPassed);

%% Rx Calibration

% Change to dB and remove the gain from the Gnu Radio.
calcPInDbShifted = 10.*log10(calcP) - rxGain;
measPInDb = calcPInDbShifted + powerShiftsForCali;

%% Antenna Nomalization
% We will deal with the antenna gain with a more dynamic model (i.e.
% compute the TX / RX LOS path gain in 3D with a constructed antenna
% pattern), instead of setting it to be a fix number.
antGain = 0;

%% Downconverter Gain
% We have instialized downConvGainInDb at the beginning of this function.

%% Final Result
pathLossInDb = txPower + antGain + downConvGainInDb - measPInDb;

end
% EOF