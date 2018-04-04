function [ polarAx, hpbw, hDediFig ] ...
    = plotAntPlanePat( polarAx, anglesInDegree, ampsInDb, FLAG_DEDI_FIG )
%PLOTANTPATTERN Plot the the antenna plane pattern in a 2D polar coordinate
%system.
%
% Inputs:
%   - polarAx
%     The polar axes to use for plotting.
%   - anglesInDegree
%     A real vector specifying the angles in degree.
%   - ampsInDb
%     The amplitude in dB.
%   - FLAG_DEDI_FIG
%     Optional. Set this to be true to generate a dedicated interactive
%     plot, instead of plotting statistically into polarAx.
%
% Yaguang Zhang, Purdue, 10/02/2017

if nargin <4
    FLAG_DEDI_FIG = false;
    hDediFig = nan;
end

try
    % If we can find the value for normalizing the antenna patter, we will
    % use it in the HPBW computation.
    maxPowerInDb = evalin('base', 'maxAntGainInDb');
    [ hpbw, anglesHpbw, powerInDbHpbw ] ...
        = computeHPBW( anglesInDegree, ampsInDb, maxPowerInDb );
catch
    [ hpbw, anglesHpbw, powerInDbHpbw ] ...
        = computeHPBW( anglesInDegree, ampsInDb );
end

if FLAG_DEDI_FIG
    % Also plot the pattern in a dedicated figure.
    hDediFig = figure;
    hold on;
    polarpattern(hDediFig, anglesInDegree, ampsInDb);
    text(0.1, 0, 'dB');
    hold off;
else
    % Plot.
    hold on;
    polarplot(deg2rad(anglesInDegree), ampsInDb-min(ampsInDb));
    polarAx.RAxis.Label.String = 'dB';
    
    polarplot(deg2rad([0 0; anglesHpbw']), ...
        [0, 0; powerInDbHpbw'-min(ampsInDb)], '--k');
    hold off;
end


end
% EOF