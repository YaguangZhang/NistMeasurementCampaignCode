function [ hFig ] = plotPdpsForOneRec(sigOutFile, F_S)
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
%
% Yaguang Zhang, Purdue, 04/30/2018

% We will only load in 1 second of the signal.
countSam = F_S; % ~ 1s of the signal.
seriesSignal = read_complex_binary(sigOutFile.name, countSam);
[~, figureSupTitle, ~] = fileparts(sigOutFile.name);

% Plot the signals. We will try to find the "tallest" bump for each
% measurement.
numPreSamples = 200;
numPostSamples = 2000;

hFig = plotOnePresentSignal(seriesSignal, ...
    numPreSamples, numPostSamples, figureSupTitle);

end
% EOF