function [ linearAmps ] = antPatDbToLinear( ampsInDb )
%ANTPATDBTOLINEAR Convert the dB amplitudes for antenna patterns to linear.
%
% In the EARS measurement campaign, we have collected the S21 amplitudes
% from both the Anzimuth and Elevation sweeps for the antenna. However, the
% resulted log files have linear amplitudes inside. This function captures
% the convertion for easy modification.
%
% Yaguang Zhang, Purdue, 10/13/2017

linearAmps = 10.^(ampsInDb./20);

end
% EOF