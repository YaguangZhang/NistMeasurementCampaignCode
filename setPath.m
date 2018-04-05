% SETPATH Add lib folders into Matlab path.
%
% Yaguang Zhang, Purdue, 03/30/2018

cd(fileparts(mfilename('fullpath')));
addpath(fullfile(pwd));
addpath(genpath(fullfile(pwd, 'lib')));

% The absolute path to the shared folder holding the data and code for the
% NIST data. Please make sure it is correct for the machine which will run
% this script.
%  - On Mac Lemma:
absPathMacLemma = '/Users/zhan1472/Google Drive/NIST Measurement Campaign';
%  - Local copy on Windows Dell:
absPathWinDell = 'C:\Users\Zyglabs\OneDrive - purdue.edu\NIST';
unknownComputerErrorMsg = ...
    ['Compute not recognized... \n', ...
    '    Please update setPath.m for your machine. '];
unknownComputerErrorId = 'setPath:computerNotKnown';
switch getenv('computername')
    case 'ZYGLABS-DELL'
        % ZYG's Dell laptop.
        ABS_PATH_TO_EARS_SHARED_FOLDER = absPathWinDell;
    case ''
        % Expected to be Lemma the Mac machine in ZYG's lab.
        assert(ismac, unknownComputerErrorMsg);
        ABS_PATH_TO_EARS_SHARED_FOLDER = absPathMacLemma;
    otherwise
        error(unknownComputerErrorId, unknownComputerErrorMsg);
end

% EOF