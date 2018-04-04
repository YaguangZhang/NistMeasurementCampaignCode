function [ ampsInDb ] = antPatLinearToDb( linearAmps )
%ANTPATLINEARTODB Convert the linear amplitudes for antenna patterns to dB.
%
% In the EARS measurement campaign, we have collected the S21 amplitudes
% from both the Anzimuth and Elevation sweeps for the antenna. However, the
% resulted log files have linear amplitudes inside. This function captures
% the convertion for easy modification.
%
% Yaguang Zhang, Purdue, 10/13/2017

ampsInDb = 20.*log10(linearAmps);

end
% EOF