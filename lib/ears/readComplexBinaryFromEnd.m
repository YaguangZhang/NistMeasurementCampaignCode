function [ samples] = readComplexBinaryFromEnd (filename, count)
% READCOMPLEXBINARYFROMEND Read in a number of samples of the complex data
% from a .out file, but from the end (instead of the start).
%
% Ref: read_complex_binary.m
%
% Yaguang Zhang, Purdue, 03/25/2018

fid = fopen (filename, 'rb');
% In a GnuRadio .out file, we have:
%     Complex - 32 bit floating point for both I and 32 Q (8 bytes in
%     total)
offset = -count.*8;
fseek(fid, offset, 'eof');
data = fread (fid, [2, count], 'float');
fclose(fid);

samples = data(1,:) + data(2,:)*1i;
[r, c] = size (samples);
samples = reshape (samples, c, r);
% EOF