function str = num2StrPrecisionForFloat(num, precision)
%NUM2STRPRECISIONFORFLOAT Convert number to string, but for float, use a
%better precision.
%
% Yaguang Zhang, Purdue, 12/27/2018

str = num2str(num);
if contains(str,'.')
    str = num2str(num, precision);   
end

end

% EOF