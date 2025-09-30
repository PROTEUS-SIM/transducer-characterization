function ScanPlane = determine_scan_plane_orientation(P,Grid,medium)
% This function computes the scan plane orientation.
%
% input
% - P: angular spectrum of hydrophone data (Nx-by-Ny-by-Nsamples)
%
% output
% - ScanPlane, struct with fields thetaX, thetaY, thetaZ
%
% See also: demo_find_angles.m
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%==========================================================================
% Parameters
%==========================================================================

% Frequency range (Hz) for finding the rotation
% One row per search range
% Search between 1.0 and 4.5 MHz, but do not include the frequencies
% between 3.2 and 3.5 MHz to filter out the spurious mode of oscillation
% (see Section III D of the arXiv preprint).
freqRange = [1.0 3.2; 3.5 4.5]*1e6;

% Estimated transducer width
W = 0.028;

% Angles of the integration vectors in k-space:
alpha = (-1:0.01:1)*pi/3;

% Search range for angle about propagation axis:
options.thetaZ = (-2:0.02:2)/180*pi;

% Don't plot the results:
options.visualize = false;

%==========================================================================
% FIND ANGLES
%==========================================================================

Fs = 1/Grid.dt;
f = (0:(Grid.Nt-1))/Grid.Nt*Fs;

% Frequency components to be included in the search:
I_search = [];
for n = 1:size(freqRange,1)
    I_search = [I_search find(f>freqRange(n,1)&f<freqRange(n,2))]; %#ok
end

A      = zeros(size(I_search)); % weights
thetaX = zeros(size(I_search));
thetaY = zeros(size(I_search));
thetaZ = zeros(size(I_search));

fprintf(' Component     | f (MHz) | X (deg) | Y (deg) | Z (deg) |\n')
fprintf('--------------------------------------------------------\n')

for n = 1:length(I_search)

    fprintf('%3d out of %3d ',n,length(I_search))
    fprintf('| %7.2f ',f(I_search(n))/1e6)
    
    % Compute the wavenumber of propagation (forward travelling waves):
    k = 2*pi*f(I_search(n))/medium.sound_speed;
    
    % Exclude integration points close to the central peak in k-space.
    % The central peak provides little information about the rotation about
    % the propagation axis.
    L = W/5;
    options.alpha = alpha(abs(alpha)>abs(2*pi/(k*L)));

    [thetaX(n), thetaY(n), thetaZ(n), A(n)] = find_angles_rectangular(...
        P(:,:,I_search(n)),k,Grid,options);

    fprintf('| %7.2f | %7.2f | %7.2f |\n',...
        thetaX(n)*180/pi,thetaY(n)*180/pi,thetaZ(n)*180/pi)

end

% Normalize the weights:
A = A/sum(A);

% Compute weighted mean of angles and convert radians to degrees:
thetaX = sum(A.*thetaX)*180/pi;
thetaY = sum(A.*thetaY)*180/pi;
thetaZ = sum(A.*thetaZ)*180/pi;

fprintf('--------------------------------------------------------\n')
fprintf(' Weighted mean |         | %7.2f | %7.2f | %7.2f |\n',...
    thetaX,thetaY,thetaZ)

ScanPlane.thetaX = thetaX;
ScanPlane.thetaY = thetaY;
ScanPlane.thetaZ = thetaZ;

end
