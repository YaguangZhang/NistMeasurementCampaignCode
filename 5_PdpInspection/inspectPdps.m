function [ ] = inspectPdps(dataTag, allSigOutFiles, absPathToSavePlot, ...
    indicesPlotsToSaveFigCopies, F_S)
%INSPECTPDPS Generate PDP .png images for signal recording .out files.
%
%   Inputs:
%       - dataTag
%         A short string to briefly characterize/describe what data are
%         beening inspected.
%       - allSigOutFiles
%         A structure array generated by rdir to specify where the .out
%         files are.
%       - absPathToSavePlot
%         The absolute path to save plots.
%       - indicesPlotsToSaveFigCopies
%         A column array with indices for the .out files specifed by
%         allSigOutFiles to be also saved as .fig copies.
%       - F_S
%         The GnuRadio sample rate for the singal recordings. The x axis
%         will be converted from sample # to time accordingly.
%
% Yaguang Zhang, Purdue, 04/30/2018

disp(['inspectPdps: Processing ', dataTag, ' data...']);
numSigOutFiles = length(allSigOutFiles);
for idxSig = 1:numSigOutFiles
    disp(['             ', num2str(idxSig), '/', num2str(numSigOutFiles)]);
    
    if exist('F_S', 'var')
        % Convert x axis unit to time if possible.
        hFig = plotPdpsForOneRec(allSigOutFiles(idxSig), F_S);
    else
        hFig = plotPdpsForOneRec(allSigOutFiles(idxSig));
        xlabel('Time (ms)');
        xticklabels(arrayfun(@(n) num2str(n/F_S*1000, '%.2f'), xticks, ...
            'UniformOutput', false));
    end
    
    plotFileName = ['PdpOverview_', dataTag, '_', num2str(idxSig)];
    
    saveas(hFig, fullfile(absPathToSavePlot, [plotFileName, '.png']));
    % Also save a .fig copy if necessary.
    if ismember(idxSig, indicesPlotsToSaveFigCopies)
        saveas(hFig, fullfile(absPathToSavePlot, [plotFileName, '.fig']));
    end

    close(hFig); 
end
%EOF