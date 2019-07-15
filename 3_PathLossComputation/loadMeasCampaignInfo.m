% LOADMEASCAMPAIGNINFO Load the measurement campaign records from
% TxInfo.txt.
%
%   When available, parameters including series number, location, TxAz,
%   TxEl, RxAz, TxEl and note, will be loaded. We will also hardcode some
%   constant parameters here.
%
%   We will use all capitalized variable names for the parameters we
%   generate here.
%
% Yaguang Zhang, Purdue, 04/09/2018

clear; clc; close all;

%% Configurations

% Add libs to current path and set ABS_PATH_TO_NIST_SHARED_FOLDER according
% to the machine name.
cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
cd('..'); setPath;

% Configure other paths accordingly.
ABS_PATH_TO_SAVE_RESULTS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'PathLossComputation');

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_RESULTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_RESULTS);
end

%% Hard Coded Parameters

% Tx Power in dBm.
TX_POWER_DBM = 23.4;
% Tx tower height in feet.
TX_HEIGHT_FEET = 5; % 60 inch.

% Rx height in meter.
RX_HEIGHT_M = 1;

% Sample rate used for GnuRadio.
F_S = 2 * 10^6;

% Signal frequency in GHz.
F_C_IN_GHZ = 28;

% Transmitter location.
TX_LAT = 39.99145872;
TX_LON = -105.2746385;

% The downconverter gain at the RX side.
DOWNCONVERTER_GAIN_IN_DB = 13.4;

% The psuedo-noise sequence generator clock frequencies in MHz for the TX
% and the RX, respectively.
TX_PN_CLK_MHz = 400;
RX_PN_CLK_MHz = 399.95;

% Set this flag to true if the PDP's are actually time-reversed. This is a
% consequence of running the sliding correlator "backwards", that is with
% the slower PN clock at the transmitter rather than the receiver. If this
% is present and true, the PDP plotting functions will inverse the input
% samples only for visalizing the PDP in the traditional way. However, path
% loss computation will not be changed as the reversed signal should have
% the same energy as the original one.
FLAG_PDP_TIME_REVERSED = true;

% Set this flag to true if it is necessary to flip the signal for the noise
% elimination procedure, too.
FLAG_PDP_TIME_REVERSED_IN_NOISE_ELI = false;

% Set this variable to be a positive value if it is necessary to get rid of
% singal samples less than a hard-coded floor.
USRP_NOISE_FLOOR_V = nan; % 0.003; % In V.

% Number of sigma to automatically set the noise floor in the noise
% elimiation process. If not set, a default value of 3 will be used.
NUM_SIGMA_FOR_THRESHOLD = 4;

%% Auto-generated Parameters

try
    TX_ALT = getElevations(TX_LAT, TX_LON);
catch
    warning(['getElevations is not able to fetch alt from Google! ', ...
        'Setting TX_ALT to the hard-coded default value...'])
    TX_ALT = 1778.9871826; % According to Google.
end

% Necessary Unit Conversion
TX_HEIGHT_M = distdim(TX_HEIGHT_FEET,'feet','meters');

SLIDE_FACTOR = TX_PN_CLK_MHz/(TX_PN_CLK_MHz-RX_PN_CLK_MHz);

%% Load Records from TxInfo.txt Files

% Find all TxInfo.txt log files.
pathToData = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, 'Data');
txInfoFilesDirs = rdir(fullfile(pathToData, '**', 'TxInfo.txt'), '', false);

% We will save the result in a column cell, corresponding to a parent dir,
% with each item being a struct array for all the measurement series under
% that folder.
[TX_INFO_LOGS, TX_INFO_LOGS_ABS_PAR_DIRS] ...
    = arrayfun(@(logDir) parseTxInfoLog(logDir.name), ...
    txInfoFilesDirs, 'UniformOutput', false);

%% Save the Results
pathFileToSaveResults = fullfile(ABS_PATH_TO_SAVE_RESULTS, ...
    'txInfoLogs.mat');

disp(' ')
disp('Saving configuration parameters...')

save(pathFileToSaveResults, 'TX_POWER_DBM', ...
    'TX_HEIGHT_FEET', 'RX_HEIGHT_M', 'TX_HEIGHT_M', ...
    'F_S', 'F_C_IN_GHZ', 'TX_LAT', 'TX_LON', 'TX_ALT', ...
    'TX_INFO_LOGS', 'TX_INFO_LOGS_ABS_PAR_DIRS', ...
    'DOWNCONVERTER_GAIN_IN_DB', 'TX_PN_CLK_MHz', 'RX_PN_CLK_MHz', ...
    'FLAG_PDP_TIME_REVERSED', 'FLAG_PDP_TIME_REVERSED_IN_NOISE_ELI', ...
    'SLIDE_FACTOR');
if exist('USRP_NOISE_FLOOR_V', 'var')
    save(pathFileToSaveResults, 'USRP_NOISE_FLOOR_V', '-append');
end
if exist('NUM_SIGMA_FOR_THRESHOLD', 'var')
    save(pathFileToSaveResults, 'NUM_SIGMA_FOR_THRESHOLD', '-append');
end

disp('    Done!')

%% Generate plotInfo.mat
% This is needed for the path loss computation scripts.

ABS_PATH_TO_DATA = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, 'Data');
ABS_PATH_TO_SAVE_PLOTS = fullfile(ABS_PATH_TO_NIST_SHARED_FOLDER, ...
    'PostProcessingResults', 'SummaryReport', 'plots');

% Create directories if necessary.
if exist(ABS_PATH_TO_SAVE_PLOTS, 'dir')~=7
    mkdir(ABS_PATH_TO_SAVE_PLOTS);
end

% Find all the parent directories for "Series_xx" data folders using regex.
disp(' ')
disp('Generating plotInfo.mat...')

% Try loading the information for the samples first.
pathToPlotInfo = fullfile(ABS_PATH_TO_SAVE_PLOTS, 'plotInfo.mat');
if exist(pathToPlotInfo, 'file')
    % The data have been processed before and the plotInfo.mat file has
    % been found. Load that to save time.
    disp('    Found plotInfo.mat; Loading the sample info in it...')
    load(pathToPlotInfo);
else
    disp('    No plotInfo.mat found; Searching for "Series" data folders...')
    % Need to actually scan the folder and find the sample folders.
    allSeriesParentDirs = rdir(fullfile(ABS_PATH_TO_DATA, '**', '*'), ...
        'regexp(name, ''(_NistFoliage$)'')');
    % Also locate all the "Series_xx" data folders for each parent
    % directory.
    allSeriesDirs = cell(length(allSeriesParentDirs),1);
    for idxPar = 1:length(allSeriesParentDirs)
        assert(allSeriesParentDirs(idxPar).isdir, ...
            ['#', num2str(idxPar), ' series parent dir should be a folder!']);
        
        curSeriesDirs = rdir( ...
            fullfile(allSeriesParentDirs(idxPar).name, '**', '*'), ...
            'regexp(name, ''(Series_\d+$)'')');
        if(isempty(curSeriesDirs))
            warning(['#', num2str(idxPar), ...
                ' series parent dir does not have any series subfolders!']);
        end
        allSeriesDirs{idxPar} = curSeriesDirs;
    end
    disp('    Saving the results...')
    % Note that the exact paths to the folders may change depending on the
    % manchine and its operating system, so only the folder names should be
    % used.
    save(pathToPlotInfo, ...
        'allSeriesParentDirs', 'allSeriesDirs');
end
disp('    Done!')

% EOF