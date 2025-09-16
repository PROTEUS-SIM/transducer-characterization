function [theta, delays] = find_angle_receive(...
    t,y,c0,Transducer,element_range,options)
%FIND_ANGLE_RECEIVE - Find the angle of the transducer with respect to an
%incoming plane wave.
%
% Returns the angle THETA in degrees and the channel delays DELAYS in
% seconds. The delays are defined such that DELAYS = 0 for the centre of
% the transducer.
%
% Input:
% - t: time array (1-by-N)
% - y: receive data (N-by-Nchannels)
% - c0: speed of sound (m/s)
% - Transducer: struct with fields:
%   - NumberOfElements
%   - Pitch
% - element_range: channels to be used for the best fit
% 
% See Fig. 6.7 of my thesis
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Find the index of the maximum of each channel
[~,idx] = max(y);

% Coordinates of the transducer elements
p = Transducer.Pitch;
N = Transducer.NumberOfElements;
x_el = (1:N)*p;
x_el = x_el - mean(x_el);

% Fit a line to the maxima
P = polyfit(x_el(element_range), t(idx(element_range)),1);

if options.showFigure
    plot((P(1)*x_el(element_range)+P(2))*1e6,element_range,...
        '-r','LineWidth',1.5)
end

% Angle of the received plane wave:
theta = -atand(P(1)*c0);

% Receive delay for each element with respect to the receive time for the
% centre of the transducer:
delays = P(1)*x_el;

end