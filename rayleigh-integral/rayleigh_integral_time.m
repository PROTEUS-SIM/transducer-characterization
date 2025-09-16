function p = rayleigh_integral_time(source, sensor, medium)
%RAYLEIGH_INTEGRAL_TIME computes the time-domain Rayleigh integral for a
%list of observation points over a set of integration points on the
%transducer surface.
%
% source, a struct with fields:
% - points:       source points [m], Npoints-by-3
% - delays:       transmit delays [s], Npoints-by-1
% - apod:         transmit apodization, real number 0 to 1, Npoints-by-1
% - samplingrate: sampling rate drive signal {Hz], real positive number
% - vdot:         transducer surface acceleration [Hz], 1-by-Ntimepoints
% - normal:       normal vector at each source point, Npoints-by-3
% - weights:      integration element surface area [m^2], Npoints-by-1
% - baffle:       'rigid', 'soft', or 'free'
%
% sensor, a struct with fields:
% - points:      sensor points [m], Npoints2-by-3
%
% medium, a struct with fields:
% - density [kg/m^3]
% - sound_speed [m/s]
% 
% RETURNS p, a real array of Npoints-by-N, where N is the number of time
% points required to capture all signal delays.
%
% See also get_transducer_integration_points, define_grid
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Developer note: This function is also used in
% microbubble-flow-simulator-gui\comparison-k-wave-analytical: Update the
% script comparison-k-wave-analytical when changing the input and output
% arguments of this function.

rho  = medium.density;
c0   = medium.sound_speed;

vdot   = source.vdot;
dA     = source.weights;
Fs     = source.samplingrate;
apod   = source.apod;
delays = source.delays;

source_points = source.points;
sensor_points = sensor.points;
normal        = source.normal;

% Compute pairwise distance matrix (Npoints2-by-Npoints1):
r1(1,:,:) = source_points; % (1-by-Npoints1-by-3 array)
r2(:,1,:) = sensor_points; % (Npoints2-by-1-by-3 array)
r = vecnorm(r1-r2,2,3);

% Obliquity factor:
beta = compute_obliquity_factor(source_points, sensor_points, ...
    normal, source.baffle);

% Time shifts in number of time samples:
shift = round((r/c0 + transpose(delays))*Fs);

% Old and new signal length:
M = length(vdot);
N = M + max(shift,[],'all');

p = zeros(size(sensor_points,1),N);

for m = 1:size(sensor_points,1)
    
    disp(['Computing pressure at observation point ' num2str(m) ' of ' ...
        num2str(size(sensor_points,1)) ' ...'])

    for n = 1:size(source_points,1)
        p(m,1+shift(m,n):M+shift(m,n)) = p(m,1+shift(m,n):M+shift(m,n))+...
            rho/(2*pi*r(m,n))*apod(n)*vdot*beta(m,n)*dA(n);
    end

end

end