% CALIBRATERX Calibartion for the received signal power.
%
% Essentially, we need to find the relationship between the "power" of the
% received valtage from Gnu Radio and the true recieved signal power.
%
% Yaguang Zhang, Purdue, 03/30/2018

clear; clc; close all;

%% Configurations

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% Configure other paths accordingly.
ABS_PATH_TO_CALI_DATA = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'Data', '20180330_Calibration');
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'Calibration');

% Set this to true if you want figures to be generated silently.
FLAG_GEN_PLOTS_SILENTLY = true;
% Set this to be true to eliminate noise (i.e. compute noise sigma) using
% signal amplitude; Otherwise, we will eliminate the noise in the real and
% imaginary parts of the signal separately.
FLAG_NOISE_ELI_VIA_AMP = true;

% Set this to be true to use the filtered (by GunRadio) outputs; Otherwise,
% the original samples without post-processing will be used.
FLAG_USE_FILTERED_OUTPUT_FILES = false;

% For the USNA measurement campaign, we have two sets of reference data,
% with the Gnu Radio gain being set to 1dB and 76dB, respectively. For the
% NIST one, we only need one set of reference data (65 dB).
rxGains = [65]; % In dB.
% Each set of data is stored in a folder named as "Gain_xxx".
calDataDirNamePrefix = 'Gain_';

% Reference received power (dBm) measured by the spectrum analyzer.
measPowers = {(-37:-2:-91)'};

% Manually ignore some of the measurements.
BOOLS_MEAS_TO_FIT = {ones(1,28)};
BOOLS_MEAS_TO_FIT{1}(end-1:end) = 0; % Ignore the last two point.

% Sample rate used for GnuRadio.
Fs = 2 * 10^6;
% For the NIST data set, we need to reverse the signal.
FLAG_PDP_TIME_REVERSED = false;

% Low pass filter for the power spectral density (PSD). 
%   - Tried before: 46000; 39500.
maxFreqPassed = 30000; % In Hz.
% High pass filter to remove the DC component.
minFreqPassed = 1; % In Hz.

% Minimum valid calculated power.
minValidCalPower = -inf; % In dB. Before: -140.

% Minimum valid estimated SNR.
minValidEstSnr = 0; % Before: 1.5.

% Number of samples to use (from the end) for each calibration recording.
% At NIST, we recorded >=8s of data for each calibration recording.
secToUsePerRecording = 3; % In seconds.
% Number of samples to discard at the beginning. We do not need to discard
% any samples at the begining for the NIST data because we've essentially
% already done so by reading in only the tails.
numStartSampsToDiscard = 0;
% After discarding these samples, furthermore only keep the middle part of
% the signal for calibration.
timeLengthAtCenterToUse = 1; % In second.

% At NIST, We only have one dataset for GunRadio gain 65 dB.
NUMS_SIGMA_FOR_THRESHOLD = [1].*4;

% For any figure generated, a .pgn screen shot is always saved; Set this to
% be true to also save a .fig version (which may dramatically slow down the
% program).
FLAG_SAVE_FIG_COPY = false;

% Parameters to overwrite when using the 60 kHz / 5 kHz LPF implemented in
% Matlab. Note: 10 kHz is optimal for the NIST implementation.
if ~FLAG_USE_FILTERED_OUTPUT_FILES
    maxFreqPassed = 10000;
end

% Regression method to use, e.g. built-in: 'robustfit', 'polyfit',
% 'regress', and external: 'linortfit2' (Orthogonal Linear Regression).
LINEAR_REGRESSION_METHOD = 'linortfit2';

%% Before Calibration

disp(' ------------- ')
disp('  calibrateRx')
disp(' ------------- ')

LogicalStr = {'false', 'true'};

disp(' ')
disp(['    FLAG_PDP_TIME_REVERSED is set to be **', ...
    LogicalStr{FLAG_PDP_TIME_REVERSED+1}, '** ...'])

% Disable the tex interpreter in figures.
set(0,'DefaultTextInterpreter','none');
% Hide figures if it's required.
if FLAG_GEN_PLOTS_SILENTLY
    set(0,'DefaultFigureVisible','off');
else
    set(0,'DefaultFigureVisible','on');
end

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
end

% Make sure the calibration data are there.
disp(' ')
disp('    Verifying calibration datasets...')

numDatasets = length(rxGains);
absPathsToCalDatasets = cell(numDatasets,1);
for idxRxGain = 1:numDatasets
    curRxGain = rxGains(idxRxGain);
    curExpFoldername = [calDataDirNamePrefix, num2str(curRxGain)];
    absPathsToCalDatasets{idxRxGain} = fullfile(ABS_PATH_TO_CALI_DATA, ...
        curExpFoldername);
    
    assert(exist( ...
        absPathsToCalDatasets{idxRxGain}, 'dir'...
        )==7, ['Expected folder ',curExpFoldername,' does not exist!']);
end

disp('    Done!')

%% Load the Calibration Datasets

disp(' ')
disp('    Loading calibration data...')

% If original .out files are used instead, we need to pass it through a 60
% kHz LPF.
if ~FLAG_USE_FILTERED_OUTPUT_FILES
    % Construct an LPF.
    Fp  = 60e3;   	% 60 kHz passband-edge frequency
    Fst = 65e3;     % Transition Width = Fst - Fp
    Ap = 0.01;      % Allowed peak-to-peak ripple
    Ast = 80;       % Stopband attenuation
    
    lpfComplex = dsp.LowpassFilter('SampleRate', Fs, ...
        'FilterType', 'FIR', 'PassbandFrequency', Fp, ...
        'StopbandFrequency', Fst, ...
        'PassbandRipple', Ap, ...
        'StopbandAttenuation', Ast ...
        );
    release(lpfComplex);
end

calData = cell(numDatasets,1);
for idxDataset = 1:numDatasets
    curDatasetAbsPath = absPathsToCalDatasets{idxDataset};
    
    if FLAG_USE_FILTERED_OUTPUT_FILES
        % Use the filtered (by GunRadio) Rx output files.
        error('The NIST data do not come with filtered recordings!');
    end
    
    % Locate the folders containing data for each flight.
    curCalDataFlightFolders = rdir(curDatasetAbsPath, ...
        'regexp(name, ''(flight_number_\d+$)'')');
    curCalDataFlightFolders = curCalDataFlightFolders(...
        [curCalDataFlightFolders.isdir]);
    
    % Make sure the log files are sorted by name before loading.
    [~, sortedIs] = sort({curCalDataFlightFolders.name});
    curCalDataFlightFolders = curCalDataFlightFolders(sortedIs);
    
    % Load all the calibration data.
    curNumMeas = length(curCalDataFlightFolders);
    curCalData = cell(curNumMeas,1);
    for idxCurMeas = 1:curNumMeas
        disp(['        Loading measurement ', num2str(idxCurMeas), ...
            '/', num2str(curNumMeas), ' for set #', num2str(idxDataset)]);
        % We will only read in the last 3 seconds of each recording.
        [curCalData{idxCurMeas}, curGlobalSettings] = ...
            parseNistFlightFolder( ...
            curCalDataFlightFolders(idxCurMeas).name, secToUsePerRecording);
        assert(curGlobalSettings.core_sample_rate == Fs, ...
            'The signal recorded for this flight is different from Fs!');
        if ~FLAG_USE_FILTERED_OUTPUT_FILES
            curCalData{idxCurMeas} = lpfComplex(curCalData{idxCurMeas});
        end
    end
    calData{idxDataset} = curCalData;
end

disp('    Done!')

%% Calibration

disp(' ')
disp('    Calibrating...')

% Compute the Rx Power (in dB) for each filtered output log file. We are
% doing the for loops again here instead of in the section above to better
% separate the code modules.
[calDataThresholded, estimatedSnrs, ... % For debugging.
    calculatedPowers] = deal(cell(numDatasets,1));
for idxDataset = 1:numDatasets
    curNumMeas = length(calData{idxDataset});
    curCalDataThr = cell(curNumMeas,1);
    curCalculatedP = nan(curNumMeas,1);
    curEstimatedSnrs = nan(curNumMeas,1);
    
    NUM_SIGMA_FOR_THRESHOLD = NUMS_SIGMA_FOR_THRESHOLD(idxDataset);
    
    for idxCurMeas = 1:curNumMeas
        % Display the progress.
        disp(['        Set: ', num2str(idxDataset),'/',...
            num2str(numDatasets), ...
            '    Measurement: ', num2str(idxCurMeas),'/',...
            num2str(curNumMeas)]);
        
        % Prepare the samples for calibration.
        curSignal = calData{idxDataset}{idxCurMeas};
        % Discard the first numStartSampsToDiscard of samples.
        curSignal = curSignal((numStartSampsToDiscard+1):end);
        % Further more, only keep the middle part for calibration.
        numSampsToKeep = ceil(timeLengthAtCenterToUse*Fs);
        numSampsCurSeries = length(curSignal);
        if numSampsToKeep > numSampsCurSeries
            warning('There are not enough samples to keep. We will use all remaining ones.');
        else
            idxRangeToKeep = floor(0.5.*numSampsCurSeries ...
                + [-1,1].*numSampsToKeep./2);
            curSignal = curSignal(idxRangeToKeep(1):idxRangeToKeep(2));
        end
        % Make sure we end up with even number of samples.
        if mod(length(curSignal),2)==1
            curSignal = curSignal(1:(end-1));
        end
        
        % Before calibration, plot the PSD of the input samples (after
        % range limitation but before noise elimination).
        X0 = curSignal;
        L0 = length(X0);
        % Apply a L-point minimum 4-term Blackman-Harris window to the
        % signal before fft.
        X0Win = blackmanharris(L0).*X0;
        % FFT results.
        Y0 = fftshift(fft(X0));
        Y0Win = fftshift(fft(X0Win));
        % Frequency domain.
        f0 = (-L0/2:L0/2-1)*(Fs/L0);
        idxDC0 = L0/2+1;
        % PSD.
        powerSpectralDen0 = abs(Y0).^2/L0;
        powerSpectralDen0Win = abs(Y0Win).^2/L0;
        
        % Plot the result in dB for debugging.
        [hPSDInputInDb, hPSDInputInDbZoomedIn] = ...
            plotPsdInDbWithLPF(powerSpectralDen0, f0, ...
            maxFreqPassed, minFreqPassed);
        set(0, 'CurrentFigure', hPSDInputInDb);
        title({['Estimated Power Spectrum Density (dB) - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)], ...
            'Before Noise Eliminiation'});
        set(0, 'CurrentFigure', hPSDInputInDbZoomedIn);
        title({['Estimated Power Spectrum Density (dB) Zoomed In - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)], ...
            'Before Noise Eliminiation'});
        [hPSDInputInDbWin, hPSDInputInDbWinZoomedIn] = ...
            plotPsdInDbWithLPF(powerSpectralDen0Win, f0, ...
            maxFreqPassed, minFreqPassed);
        set(0, 'CurrentFigure', hPSDInputInDbWin);
        title({['Estimated Power Spectrum Density (dB) - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)], ...
            'Before Noise Eliminiation, Windowed'});
        set(0, 'CurrentFigure', hPSDInputInDbWinZoomedIn);
        title({['Estimated Power Spectrum Density (dB) Zoomed In - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)], ...
            'Before Noise Eliminiation, Windowed'});
        
        % For calibration, first of all, threshold the I&Q waveforms of
        % each signal vector to eliminate corss-correlation and system
        % noise.
        [signalReal, ~, hNoiseSigmaReal] = ...
            thresholdWaveform(real(curSignal), true);
        [signalImag, ~, hNoiseSigmaImag] = ...
            thresholdWaveform(imag(curSignal), true);
        % Update: it makes more sense to eliminate pts according to the
        % amplitude of the signal, instead of doing it separately for the
        % real and image parts.
        [~, boolsEliminatedPts, hNoiseSigmaAmp] = ...
            thresholdWaveform(abs(curSignal), true);
        curSignalEliminated = curSignal;
        curSignalEliminated(boolsEliminatedPts) = 0;
        
        % Depending on the flag, choose to use which version of
        % noise-eliminated signal (via real and imaginary parts separately,
        % or via the amplitude as a whole).
        if FLAG_NOISE_ELI_VIA_AMP
            curCalDataThr{idxCurMeas} = curSignalEliminated;
        else
            % Compute the complex FFT of the resulted signal; Store it in
            % the cell curCalDataThr for debugging.
            curCalDataThr{idxCurMeas} = signalReal+1i.*signalImag;
        end
        
        % Also get rid of everything below the USRP noise floor if
        % USRP_NOISE_FLOOR_V is specified in the base workspace.
        if evalin('base','exist(''USRP_NOISE_FLOOR_V'', ''var'')')
            USRP_NOISE_FLOOR_V = evalin('base', 'USRP_NOISE_FLOOR_V');
            curCalDataThr{idxCurMeas}...
                (abs(curCalDataThr{idxCurMeas})<USRP_NOISE_FLOOR_V) = 0;
        end

        % Signal (noise eliminiated) to process.
        X = curCalDataThr{idxCurMeas};
        L = length(X);
        
        % Apply a L-point minimum 4-term Blackman-Harris window to the
        % signal before fft.
        %  X = blackmanharris(L).*X;
        % FFT results.
        Y = fftshift(fft(X));
        % Frequency domain.
        f = ((-L/2:L/2-1)*(Fs/L))';
        idxDC = L/2+1;
        % PSD.
        powerSpectralDen = abs(Y).^2/L;
        
        % Plot the result for debugging.
        hSigOverview = figure; hold on;
        numSamplesToPlot = 2000000; % For 2 MSamples/s signal => 1 s.
        plot(1:numSamplesToPlot, abs(curSignal(1:numSamplesToPlot)));
        axis tight; grid on;
        title('Signal Overview');
        xlabel('Sample #'); ylabel('Amplitude');         
        
        hPSD = figure; hold on;
        hPowerSpectralDen = plot(f,powerSpectralDen);
        curAxis = axis; curAxis(1) = f(1); curAxis(2) = f(end);
        if ~isinf(maxFreqPassed)
            hLPF = plot([maxFreqPassed, -maxFreqPassed; ...
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
            hHPF = plot([minFreqPassed, -minFreqPassed; ...
                minFreqPassed, -minFreqPassed], ...
                [curAxis(3),curAxis(3); ...
                curAxis(4),curAxis(4)], '-.k');
            x = [-minFreqPassed minFreqPassed ...
                minFreqPassed -minFreqPassed];
            y = [curAxis(3) curAxis(4) curAxis(4) curAxis(3)];
            patch(x,y,[1,1,1].*0.6,'FaceAlpha',0.3,'LineStyle','none');
        end
        
        % Compute the power.
        boolsFPassed = abs(f)<=maxFreqPassed ...
            & abs(f)>=minFreqPassed;
        % Compute the power by integral. Note that we will always discard
        % the DC component here (although it may be passed by the filters).
        psdPassed = powerSpectralDen;
        psdPassed(~boolsFPassed) = 0;
        psdPassed(idxDC) = 0;
        curCalculatedP(idxCurMeas) = trapz(f, psdPassed);
        % For the noise, we compute the power right outside of the LPF but
        % limit the integral range to be as wide as the LPF. That is, as
        % the reference noise power, compute the power from maxFreqPassed
        % to 2*maxFreqPassed in both the positive and negative parts.
        boolsFFiltered = abs(f)<=2*maxFreqPassed & (~boolsFPassed);
        psdFiltered = powerSpectralDen;
        psdFiltered(~boolsFFiltered) = 0;
        psdFiltered(idxDC) = 0;
        powerFiltered = trapz(f, psdFiltered);
        curEstimatedSnrs(idxCurMeas) = ...
            curCalculatedP(idxCurMeas)/powerFiltered;
        
        text(min(maxFreqPassed, 50000), mean([curAxis(3) curAxis(4)]), ...
            ['Estimated SNR = ', ...
            num2str(curEstimatedSnrs(idxCurMeas), '%.2f')]);
        hold off;
        title(['Estimated Power Spectrum Density - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)]);
        xlabel('f (Hz)'); ylabel('Estimated PSD (V^2/Hz)'); axis(curAxis);
        if ~isinf(maxFreqPassed) && minFreqPassed>0
            legend([hPowerSpectralDen, hLPF(1), hHPF(1)], ...
                'PSD', 'LPF', 'HPF');
        elseif isinf(maxFreqPassed) && minFreqPassed>0
            legend([hPowerSpectralDen, hHPF(1)], ...
                'PSD', 'HPF');
        elseif ~isinf(maxFreqPassed) && minFreqPassed<=0
            legend([hPowerSpectralDen, hLPF(1)], ...
                'PSD', 'LPF');
        else
            legend(hPowerSpectralDen, 'PSD');
        end
        transparentizeCurLegends;
        grid minor;
        
        % The same PSD plot in dB.
        [hPSDdB, hPSDdBZoomedIn] = plotPsdInDbWithLPF(powerSpectralDen, f, maxFreqPassed, ...
            minFreqPassed, curEstimatedSnrs(idxCurMeas));
        set(0, 'CurrentFigure', hPSDdB);
        title(['Estimated Power Spectrum Density (dB) - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)]);
        set(0, 'CurrentFigure', hPSDdBZoomedIn);
        title(['Estimated Power Spectrum Density (dB) Zoomed In - Set #', ...
            num2str(idxDataset), ' Pt #', num2str(idxCurMeas)]);
        
        % Paths to save the plots.
        pathSigOverviewToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), ...
            '-sig-overview']);
        pathInputPSDdBFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), ...
            '-psd-input-dB']);
        pathWindowedInputPSDdBFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), ...
            '-psd-input-windowed-dB']);
        pathNewPSDFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), '-psd']);
        pathNewPSDdBFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), '-psd-dB']);
        pathNewCompNoiseSigmaToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
            ['set-',num2str(idxDataset),'-pt-',num2str(idxCurMeas), '-noise-sigma-']);
        % A .png figure for easy access.
        saveas(hSigOverview, [pathSigOverviewToSave, '.png']);        
        saveas(hPSDInputInDb, [pathInputPSDdBFileToSave, '.png']);
        saveas(hPSDInputInDbZoomedIn, [pathInputPSDdBFileToSave, '-zoomed-in.png']);
        saveas(hPSDInputInDbWin, [pathWindowedInputPSDdBFileToSave, ...
            '.png']);
        saveas(hPSDInputInDbWinZoomedIn, [pathWindowedInputPSDdBFileToSave, ...
            '-zoomed-in.png']);
        saveas(hPSD, [pathNewPSDFileToSave, '.png']);
        saveas(hPSDdB, [pathNewPSDdBFileToSave, '.png']);
        saveas(hPSDdBZoomedIn, [pathNewPSDdBFileToSave, '-zoomed-in.png']);
        % Also a .fig copy.
        if FLAG_SAVE_FIG_COPY
            saveas(hSigOverview, [pathSigOverviewToSave, '.fig']);   
            saveas(hPSDInputInDb, [pathInputPSDdBFileToSave, '.fig']);
            saveas(hPSDInputInDbZoomedIn, [pathInputPSDdBFileToSave, '-zoomed-in.fig']);
            saveas(hPSDInputInDbWin, [pathWindowedInputPSDdBFileToSave, ...
                '.fig']);
            saveas(hPSDInputInDbWinZoomedIn, [pathWindowedInputPSDdBFileToSave, ...
                '-zoomed-in.fig']);
            saveas(hPSD, [pathNewPSDFileToSave, '.fig']);
            saveas(hPSDdB, [pathNewPSDdBFileToSave, '.fig']);
            saveas(hPSDdBZoomedIn, [pathNewPSDdBFileToSave, '-zoomed-in.fig']);
        end
        
        % Close the figures if we do not need to show them.
        if FLAG_GEN_PLOTS_SILENTLY
            close([hSigOverview hPSDInputInDb hPSDInputInDbZoomedIn ...
                hPSDInputInDbWin hPSDInputInDbWinZoomedIn ...
                hPSD hPSDdB hPSDdBZoomedIn]);
        end
        
        % Plot the noise elimination for the real & imaginary parts only if
        % we use them for calculating the signal power.
        if ~FLAG_NOISE_ELI_VIA_AMP
            saveas(hNoiseSigmaReal, [pathNewCompNoiseSigmaToSave, 'real.png']);
            saveas(hNoiseSigmaImag, [pathNewCompNoiseSigmaToSave, 'imag.png']);
            if FLAG_SAVE_FIG_COPY
                saveas(hNoiseSigmaReal, [pathNewCompNoiseSigmaToSave, 'real.fig']);
                saveas(hNoiseSigmaImag, [pathNewCompNoiseSigmaToSave, 'imag.fig']);
            end
            
            % Close the figures if we do not need to show them.
            if FLAG_GEN_PLOTS_SILENTLY
                close([hNoiseSigmaReal hNoiseSigmaImag]);
            end
        end
        
        % Always plot the amplitude noise elimination for comparsion.
        saveas(hNoiseSigmaAmp, [pathNewCompNoiseSigmaToSave, 'amp.png']);
        if FLAG_SAVE_FIG_COPY
            saveas(hNoiseSigmaAmp, [pathNewCompNoiseSigmaToSave, 'amp.fig']);
        end
        
        % Close the figures if we do not need to show them.
        if FLAG_GEN_PLOTS_SILENTLY
            close(hNoiseSigmaAmp);
        end
    end
    calDataThresholded{idxDataset} = curCalDataThr;
    % Change to dB and remove the gain from the Gnu Radio.
    calculatedPowers{idxDataset} = 10.*log10(curCalculatedP) ...
        - rxGains(idxDataset);
    estimatedSnrs{idxDataset} = curEstimatedSnrs;
end

disp('    Done!')

%% Plot

disp(' ')
disp('    Plotting...')

seriesColors = parula;
[numSeriesColors, ~] = size(seriesColors);
rng(2);
indicesColorToUse = randi([1 numSeriesColors],1,numDatasets);

% Calculated power vs measured power.
hFigCalibrationCalcVsMeas = figure; hold on;
calPs = vertcat(calculatedPowers{:});
calPs = calPs(~isinf(calPs));
meaPs = vertcat(measPowers{:});
axisToSet = [min(meaPs) max(meaPs) ...
    min(calPs) max(calPs)];
% Fitted least-squares line.
[hLsLines, hLsInvLines, lsLinesPolys, lsLinesPolysInv] = ...
    deal(cell(numDatasets,1));
for idxDataset = 1:numDatasets
    xs = measPowers{idxDataset};
    ys = calculatedPowers{idxDataset};
    
    % For plotting, only avoid points with inf as calculated power.
    xsToShow = xs(~isinf(ys));
    ysToShow = ys(~isinf(ys));
    
    % Plot the results.
    colorToUse = seriesColors(indicesColorToUse(idxDataset),:);
    % Non-inf points.
    scatter(xsToShow, ysToShow, '*', 'MarkerEdgeColor', colorToUse, ...
        'LineWidth',1.5);
end
% Set the visible area of the plot now according to the data points shown.
axis(axisToSet); axis equal; finalAxis = axis; axis manual;
% Add the lslines.
[fittedMeaPs, fittedCalPs] = deal(cell(numDatasets,1));
for idxDataset = 1:numDatasets
    xs = measPowers{idxDataset};
    ys = calculatedPowers{idxDataset};
    
    % For fitting, remove points with too low estimated SNR.
    boolsPtsToFit = estimatedSnrs{idxDataset}>=minValidEstSnr;
    % Also get rid of measurements to use in line fitting.
    boolsPtsToFit = boolsPtsToFit & BOOLS_MEAS_TO_FIT{idxDataset}';
    xsToFit = xs(boolsPtsToFit);
    ysToFit = ys(boolsPtsToFit);
    
    % For fitting, remove points with too low calculated power.
    xsToFit = xsToFit(xsToFit>=minValidCalPower);
    ysToFit = ysToFit(ysToFit>=minValidCalPower);
    
    % Record the data points used in the fitting process.
    fittedMeaPs{idxDataset} = xsToFit;
    fittedCalPs{idxDataset} = ysToFit;
    
    % Cover the unused points in the plot.
    hIgnoredPts = plot(xs(~boolsPtsToFit),ys(~boolsPtsToFit), ...
        'r*', 'LineWidth',1.5);
    
    switch LINEAR_REGRESSION_METHOD
        case 'polyfit'
            % Linear fitting.
            lsLinePoly = polyfit(xsToFit, ysToFit, 1);
            % For future use & comparison.
            lsLinePolyInv = polyfit(ysToFit, xsToFit, 1);
        case 'robustfit'
            lsLinePoly = robustfit(xsToFit, ysToFit);
            lsLinePolyInv = robustfit(ysToFit, xsToFit);
            % To go with the output order for polyfit.
            lsLinePoly = lsLinePoly(end:-1:1);
            lsLinePolyInv = lsLinePolyInv(end:-1:1);
        case 'regress'
            lsLinePoly = regress(ysToFit,[ones(length(xsToFit),1) xsToFit]);
            lsLinePolyInv = regress(xsToFit,[ones(length(ysToFit),1) ysToFit]);
            % To go with the output order for polyfit.
            lsLinePoly = lsLinePoly(end:-1:1);
            lsLinePolyInv = lsLinePolyInv(end:-1:1);
        case 'linortfit2'
            % Linear fitting.
            lsLinePoly = linortfit2(xsToFit, ysToFit);
            % For future use & comparison.
            lsLinePolyInv = linortfit2(ysToFit, xsToFit);
        otherwise
            error('Linear regression method not supported!')
    end
    
    xRangeToShow = linspace(finalAxis(1),finalAxis(2));
    % First plot the cooresponding lsLinePolyInv as background, i.e. we
    % will plot inverse(lsLinePolyInv) here.
    hLsInvLines{idxDataset} = plot(xRangeToShow, ...
        polyval(...
        [1/lsLinePolyInv(1), -lsLinePolyInv(2)/lsLinePolyInv(1)], ...
        xRangeToShow), ...
        'LineStyle', '-.', 'Color',ones(1,3).*0.8);
    % Then plot the fitted line.
    valuesLsLine = polyval(lsLinePoly,xRangeToShow);
    hLsLines{idxDataset} = ...
        plot(xRangeToShow,valuesLsLine, ...
        'Color',colorToUse,'LineStyle', '--');
    
    % Show the polynomial on the plot.
    if lsLinePoly(2)>0
        strPoly=['y = ',num2str(lsLinePoly(1)),'x+',num2str(lsLinePoly(2))];
    elseif lsLinePoly(2)<0
        strPoly=['y = ',num2str(lsLinePoly(1)),'x',num2str(lsLinePoly(2))];
    else % lsLinePoly(2)==0
        strPoly=['y = ',num2str(lsLinePoly(1)),'x'];
    end
    idxMiddlePtToFit = floor(length(xsToFit)/2);
    % Black bold copy as background for clarity.
    text(xsToFit(idxMiddlePtToFit), ...
        ysToFit(idxMiddlePtToFit), strPoly, ...
        'Rotation', rad2deg(atan(lsLinePolyInv(1))), ...
        'FontWeight', 'bold', 'Color', 'white', ...
        'VerticalAlignment', 'top');
    text(xsToFit(idxMiddlePtToFit), ...
        ysToFit(idxMiddlePtToFit), strPoly, ...
        'Rotation',rad2deg(atan(lsLinePolyInv(1))), ...
        'Color', seriesColors(indicesColorToUse(idxDataset),:), ...
        'VerticalAlignment', 'top');
    
    lsLinesPolys{idxDataset} = lsLinePoly;
    lsLinesPolysInv{idxDataset} = lsLinePolyInv;
end
if ~isinf(minValidCalPower)
    plot([finalAxis(1) finalAxis(2)], [minValidCalPower minValidCalPower], ...
        'r--');
    x = [finalAxis(1) finalAxis(1) finalAxis(2) finalAxis(2)];
    minY = -200;
    y = [minY minValidCalPower minValidCalPower minY];
    patch(x,y,[1,1,1].*0.6,'FaceAlpha',0.3,'LineStyle','none');
end
title('Calibration results');
xlabel('Measured Power (dB)');
ylabel('Calculated Power (dB)');
grid minor; hold off;

% Similarly, measured power vs calculated power.
hFigCalibrationMeasVsCalc = figure; hold on;
axisToSet = [ min(calPs) max(calPs) ...
    min(meaPs) max(meaPs)];
% Calibration points.
for idxDataset = 1:numDatasets
    xs = calculatedPowers{idxDataset};
    ys = measPowers{idxDataset};
    
    % % For plotting, only show points to fit.
    %  boolsPtsToFit = BOOLS_MEAS_TO_FIT{idxDataset}';
    % xsToShow = xs((~isinf(ys))&boolsPtsToFit);
    %  ysToShow = ys((~isinf(ys))&boolsPtsToFit);
    
    % For plotting, only avoid points with inf as calculated power.
    xsToShow = xs(~isinf(ys));
    ysToShow = ys(~isinf(ys));
    
    % Plot the results.
    colorToUse = seriesColors(indicesColorToUse(idxDataset),:);
    % Non-inf points.
    scatter(xsToShow, ysToShow, '*', 'MarkerEdgeColor', colorToUse, ...
        'LineWidth',1.5);
end
% Set the visible area of the plot now according to the data points shown.
axis(axisToSet); axis equal; finalAxis = axis; axis manual;
% Add the lslines.
for idxDataset = 1:numDatasets
    xs = calculatedPowers{idxDataset};
    ys = measPowers{idxDataset};
    
    % For fitting, remove points with too low estimated SNR.
    boolsPtsToFit = estimatedSnrs{idxDataset}>=minValidEstSnr;
    % Also get rid of measurements to use in line fitting.
    boolsPtsToFit = boolsPtsToFit & BOOLS_MEAS_TO_FIT{idxDataset}';
    xsToFit = xs(boolsPtsToFit);
    ysToFit = ys(boolsPtsToFit);
    
    % For fitting, remove points with too low calculated power.
    xsToFit = xsToFit(xsToFit>=minValidCalPower);
    ysToFit = ysToFit(ysToFit>=minValidCalPower);
    
    % Cover the unused points in the plot.
    hIgnoredPts = plot(xs(~boolsPtsToFit),ys(~boolsPtsToFit), ...
        'r*', 'LineWidth',1.5);
    
    xRangeToShow = linspace(finalAxis(1),finalAxis(2));
    % Plot the fitted line.
    lsLinePolyInv = lsLinesPolysInv{idxDataset};
    valuesLsLine = polyval(lsLinePolyInv, xRangeToShow);
    hLsLines{idxDataset} = ...
        plot(xRangeToShow,valuesLsLine, ...
        'Color',colorToUse,'LineStyle', '--');
    
    % Show the polynomial on the plot.
    if lsLinePolyInv(2)>0
        strPoly=['y = ',num2str(lsLinePolyInv(1)),'x+',num2str(lsLinePolyInv(2))];
    elseif lsLinePolyInv(2)<0
        strPoly=['y = ',num2str(lsLinePolyInv(1)),'x',num2str(lsLinePolyInv(2))];
    else % lsLinePolyInv(2)==0
        strPoly=['y = ',num2str(lsLinePolyInv(1)),'x'];
    end
    idxMiddlePtToFit = floor(length(xsToFit)/2);
    % Black bold copy as background for clarity.
    text(xsToFit(idxMiddlePtToFit), ...
        ysToFit(idxMiddlePtToFit), strPoly, ...
        'Rotation', rad2deg(atan(lsLinePolyInv(1))), ...
        'FontWeight', 'bold', 'Color', 'white', ...
        'VerticalAlignment', 'top');
    text(xsToFit(idxMiddlePtToFit), ...
        ysToFit(idxMiddlePtToFit), strPoly, ...
        'Rotation',rad2deg(atan(lsLinePolyInv(1))), ...
        'Color', seriesColors(indicesColorToUse(idxDataset),:), ...
        'VerticalAlignment', 'top');
end
if ~isinf(minValidCalPower)
    plot([finalAxis(1) finalAxis(2)], [minValidCalPower minValidCalPower], ...
        'r--');
    x = [finalAxis(1) finalAxis(1) finalAxis(2) finalAxis(2)];
    minY = -200;
    y = [minY minValidCalPower minValidCalPower minY];
    patch(x,y,[1,1,1].*0.6,'FaceAlpha',0.3,'LineStyle','none');
end
title('Calibration results');
xlabel('Calculated Power (dB)');
ylabel('Measured Power (dB)');
grid minor; hold off;

% Save the figures. We will embed key parameters into the figure file name,
% e.g.:
%   Calibration_20170915_center_1s_matlabLPF_20kHz_ignore_2_1st_1_2nd.png
if FLAG_USE_FILTERED_OUTPUT_FILES
    matlabLPFStr = '';
else
    matlabLPFStr = ['matlabLPF_', num2str(Fp/1000), 'kHz_'];
end
[ignoreStr, ignoreStrSet1, ignoreStrSet2] = deal('');
% For the NIST dataset, we only have one calibration line.
if ~all([BOOLS_MEAS_TO_FIT{1:end}])
    ignoreStr = '_ignore_';
    if ~all(BOOLS_MEAS_TO_FIT{1})
        ignoreStrSet1 = [num2str(sum(~BOOLS_MEAS_TO_FIT{1})), '_1st'];
    end
    % if ~all(BOOLS_MEAS_TO_FIT{2})
    %     ignoreStrSet2 = ['_', num2str(sum(~BOOLS_MEAS_TO_FIT{2})), '_2nd'];
    % end
end
pathCalFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, ...
    ['Calibration_', datestr(datetime('now'), 'yyyymmdd'), '_', ...
    LINEAR_REGRESSION_METHOD, ...
    '_ths_', num2str(NUMS_SIGMA_FOR_THRESHOLD(1)), ... % '_', num2str(NUMS_SIGMA_FOR_THRESHOLD(2)), ...
    '_center_', num2str(timeLengthAtCenterToUse), 's_', ...
    matlabLPFStr, ...
    'range_', ...
    num2str(minFreqPassed), 'Hz', ...
    num2str(maxFreqPassed/1000), 'kHz_', ...
    ignoreStr, ignoreStrSet1, ignoreStrSet2]);
saveas(hFigCalibrationCalcVsMeas, [pathCalFileToSave, '_CalcVsMeas.png']);
saveas(hFigCalibrationMeasVsCalc, [pathCalFileToSave, '_MeasVsCal.png']);
if FLAG_SAVE_FIG_COPY
    saveas(hFigCalibrationCalcVsMeas, [pathCalFileToSave, '_CalcVsMeas.fig']);
    saveas(hFigCalibrationMeasVsCalc, [pathCalFileToSave, '_MeasVsCal.fig']);
end

if FLAG_GEN_PLOTS_SILENTLY
    set(0,'DefaultFigureVisible','on');
end

% Save the calibration points and the resulted polynomials for the fitted
% lines.
pathCalFileToSave = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'lsLinesPolys');
save([pathCalFileToSave, '.mat'], ...
    'lsLinesPolys', 'lsLinesPolysInv', 'fittedMeaPs', 'fittedCalPs', ...
    'rxGains');

if FLAG_GEN_PLOTS_SILENTLY
    close all;
end

disp('    Done!')

% EOF