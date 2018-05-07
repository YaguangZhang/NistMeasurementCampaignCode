function [ samples] = readComplexBinaryInRange (filename, range)
%COUNTCOMPLEXBINARY Read only the complex sample points in a GnuRadio .out
%file in the range of (range(1):range(2)).
%
% Ref: read_complex_binary.m
%
% Yaguang Zhang, Purdue, 10/06/2017

fid = fopen (filename, 'rb');
% In a GnuRadio .out file, we have:
%     Complex - 32 bit floating point for both I and 32 Q (8 bytes in
%     total per sample).
offset = (range(1)-1).*8;
fseek(fid, offset, 'bof');
data = fread (fid, [2, range(2)-range(1)+1], 'float');
fclose(fid);

samples = data(1,:) + data(2,:)*1i;
[r, c] = size (samples);
samples = reshape (samples, c, r);

end
% EOF