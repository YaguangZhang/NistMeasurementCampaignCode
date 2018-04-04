function [ numComplexSamples ] = countComplexBinary (filename)
%COUNTCOMPLEXBINARY Count the number of complex sample points in a GnuRadio
%.out file.
%
% Ref: read_complex_binary.m
%
% Yaguang Zhang, Purdue, 10/06/2017

fid = fopen (filename, 'rb');
fseek(fid, 0, 'eof');
position = ftell(fid);    % In byte.
fclose(fid);

% In a GnuRadio .out file, we have:
%     Complex - 32 bit floating point for both I and 32 Q

numComplexSamples = position/8;

end
% EOF