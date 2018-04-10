function [ amps ] = antPatInter(patAz, patEl, azs, els, INTER_METHOD)
%ANTPATINTER Antenna pattern interpolation.
%
% Inputs:
%   - patAz, patEl
%     The antenna patterns, for the Azimuth and Elevation sweeps,
%     respectively; Each of which is a struct containing fields:
%       - azs
%         The azimuth angles in degree from set [0, 360).
%       - els
%         The elevation angles in degree from set [0, 360).
%       - amps
%         The linear amplitudes of the samples.
%       - phases
%         The phases of the samples.
%     All of these fields contains a column vector with each row
%     corresponding to a sweep sample.
%   - azs, els
%     The points that needed to interpolate. The ranges of them are not
%     limited, but the unit for them should be degree.
% Output:
%   - amps
%     Interpolated amplitude (linear).
%
% For now, we simply apply a customized linear interpolation.
%
% Yaguang Zhang, Purdue, 10/02/2017

% Methods supported: 'OnLine' and 'WeightedSum'.
if nargin < 5
    INTER_METHOD = 'OnLine';
end

assert(all(size(azs) == size(els)), ...
    'The sizes of inputs azs and els should be the same!');

switch INTER_METHOD
    case 'OnLine'
        %% Linearly Interpolate on the (Azimuth, Elevation) Plane
        %
        % Good! Converges to the measurements.
        %
        % BAD! Too ad-hoc! Shape is way off expectation. This method uses
        % reference points that may be far away from the input direction.
        % Also, amplitude result may change for the same direction with
        % different forms, e.g. (30, 30) vs. (210, 150) in terms of (az,
        % el).
        
        % Limit the ranges of input angles to [0, 360).
        azs = mod(azs, 360);
        els = mod(els, 360);
        
        % Draw a line with slope -1 through the input (az, el) to find the
        % cooresponding crossed points on the x(az) and y(el) axes, i.e.
        %    el = -az + b
        % Note that we only need to worry about the angles ranging from 0
        % to 360 degrees for azs and els, but bs will be ranging from 0 to
        % 720.
        bs = azs+els;
        
        % Fetch the amplitude for the points on the axis via linear
        % interpolation.
        amp0sAz = interp1([patAz.azs; patAz.azs(2:end)+360], ...
            [patAz.amps; patAz.amps(2:end)], bs);
        amp0sEl = interp1([patEl.els; patEl.els(2:end)+360], ...
            [patEl.amps; patEl.amps(2:end)], bs);
        
        % Find all the 0 b's.
        boolsToSkip = bs(:) == 0;
        indicesToFit = 1:numel(azs);
        indicesToFit = indicesToFit(~boolsToSkip);
        
        % Now linearly interpolate, along the line we drew, between the
        % fetched value, i.e.
        %    (EL0 to P) / (EL0 to AZ0) = az / b
        %  = (amp - ampEl0) / (ampAz0 ampEl0)
        amps = nan(size(azs));
        amps(indicesToFit) = arrayfun(@(idx) interp1( [0, bs(idx)], ...
            [amp0sEl(idx), amp0sAz(idx)], ...
            azs(idx)), indicesToFit);
        
        % For zero b's, we need to return amp0 for AZ = EL = 0. We will
        % just use the value for zero az from patAz.
        assert(patAz.azs(1) == 0, 'Samples in patAz should start with az=0!');
        amps(boolsToSkip) = patAz.amps(1);
    case 'WeightedSum'
        %% Weighted Sum Method
        %
        % Good! Fallows the expected shape better. Also converges to the
        % measurements.
        
        % We will take advantage here the Matlab convention about Azimuth
        % and Elevation:
        %
        %   The azimuth angle of a vector is the angle between the x-axis
        %   and the orthogonal projection of the vector onto the xy plane.
        %   The angle is positive in going from the x axis toward the y
        %   axis. Azimuth angles lie between –180 and 180 degrees. The
        %   elevation angle is the angle between the vector and its
        %   orthogonal projection onto the xy-plane. The angle is positive
        %   when going toward the positive z-axis from the xy plane. These
        %   definitions assume the boresight direction is the positive
        %   x-axis.
        %
        
        % Project the directions we want to x-y and x-z planes,
        % respectively, to get the angles for fetching the antenna pattern
        % data as reference. Note that the angle we need for the x-y plane
        % is azs, which is already available, but it may not be of the
        % correct value corresponding to els0.
        [xs, ys, zs] = sph2cart(deg2rad(azs), deg2rad(els), ones(size(azs)));
        rCut = (1-xs.^2).^0.5;
        azs0 = mod(rad2deg(atan2(sign(ys).*rCut, xs)), 360);
        els0 = mod(rad2deg(atan2(sign(zs).*rCut, xs)), 360);
        
        % Fetch the corresponding amplitude for (0, el) and (az, 0).
        amp0sAz = interp1(patAz.azs, patAz.amps, azs0);
        amp0sEl = interp1(patEl.els, patEl.amps, els0);
        % Compute how close the point is to plane x-y and x-z, and use the
        % results to get weights for averaging two reference amplitudes,
        % amp2sAz and amp0sEl. Note:
        %     -  360 degree is essentially 0, too;
        %      - We have weights from [0, 1];
        %     -  And the closer the point is towards the x-y / x-z plane,
        %     the smaller abs(z) / abs(y) will become, and the smaller the
        %     weight will be for amp0sEl / amp2sAz.
        wAzs = abs(ys);
        wEls = abs(zs);
        amps = (amp0sAz .* wAzs + amp0sEl .* wEls)./(wAzs + wEls);
        % For zero b's, we need to return amp0 for AZ = EL = 0. We will
        % just use the value for zero az from patAz.
        assert(patAz.azs(1) == 0, 'Samples in patAz should start with az=0!');
        amps(azs0==0&els0==0) = patAz.amps(1);
        assert(all(~isnan(amps(:))), ...
            'Resulted amps should only contain non-nan elements!');
        
    otherwise
        error(['Unsupported antenna pattern interpolation method: ', INTER_METHOD, '!'])
        
end

end
% EOF