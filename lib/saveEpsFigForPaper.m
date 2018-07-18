function [ ] = saveEpsFigForPaper(hFig, fullPathToSaveFig)
%SAVEEPSFIGFORPAPER Save an .eps version for figure hFig at
%fullPathToSaveFig.
%
% Note: In order to use export_fig, one needs to install ghostscript at
%       http://www.ghostscript.com 
% and manually locate pdftops (if asked) at 
%       NistMeasurementCampaignCode/lib/ext/xpdf-tools-win-4.00
%
% Yaguang Zhang, Purdue, 07/16/2018

[dirToSave, figName, fileExt] = fileparts(fullPathToSaveFig);
% Create directories if necessary.
if exist(dirToSave, 'dir')~=7
    mkdir(dirToSave);
end
if ~strcmpi(fileExt, '.eps')
    warning( ...
        'The file extention specified is not .eps and will be ignored.')
end

epsFullPathToSave = fullfile(dirToSave, [figName, '.eps']);

curFigure = gcf;

set(0, 'currentfigure', hFig);
export_fig(epsFullPathToSave, '-eps', '-transparent'); 

set(0, 'currentfigure', curFigure);
end

% EOF