function [thetaX, thetaY, thetaZ, Imax] = find_angles_rectangular(...
    A,k,Grid,options)
%FIND_ANGLES_RECTANGULAR finds the rotation angles that define the rotation
%of the sensor plane with respect to the source plane assuming a
%rectangular aperture.
%
% A is a matrix that represent the angular spectrum of the sensor data at a
% single frequency corresponding to wavenumber k.
%
% Grid should have the following fields: Nx, Ny, dx, and dy.
% The z-axis of Grid should corresponds to the axis of propagation.
%
% Returns the angles thetaX, thetaY, and thetaZ in radians.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Enlarge spatial grid, resulting in a more finely spaced k-space grid,
% allowing for a more accurate angle search:
expansion_factor = 20;
Nx = Grid.Nx*expansion_factor;
Ny = Grid.Ny*expansion_factor;

% Set up k-space grid:
kx  = 2*pi*(ceil(-Nx/2):ceil(Nx/2-1))/(Nx*Grid.dx);
ky  = 2*pi*(ceil(-Ny/2):ceil(Ny/2-1))/(Ny*Grid.dy);

% Interpolate the data:
A = sinc_interpolation_periodic(Grid.kx,A,kx);
A = transpose(A);
A = sinc_interpolation_periodic(Grid.ky,A,ky);
A = transpose(A);
I = abs(A.^2);

Imax = max(I(:));

%==========================================================================
% FIND ANGLES
%==========================================================================

% Define integration points in k-space. These are points in k-space where
% the angular spectrum is expected to be strong. For a rectangular source
% aperture, the angular spectrum will be strong for points along the kx and
% ky axes.
if isfield(options,'alpha')
    alpha = options.alpha;
else
    alpha = (-1:0.02:1)*pi/8;
end
v1 = [k*sin(alpha); 0*sin(alpha); k*cos(alpha)];
v2 = [0*sin(alpha); k*sin(alpha); k*cos(alpha)];
v  = [v1'; v2'];
kxi = v(:,1);
kyi = v(:,2);

% Search range for the angle about the propagation axis:
if isfield(options,'thetaZ')
    thetaZ = options.thetaZ;
else
    thetaZ = (-15:0.01:15)/180*pi;
end

[thetaX,thetaY,thetaZ] = find_angles(I,thetaZ,kx,ky,kxi,kyi,k);

%==========================================================================
% SHOW ROTATED INTEGRATION POINTS
%==========================================================================

if isfield(options,'visualize') && ~options.visualize
    return
end

% Get rotation matrices:
Rx = get_rotation_matrix(thetaX,1);
Ry = get_rotation_matrix(thetaY,2);
Rz = get_rotation_matrix(thetaZ,3);

R = Rz*Ry*Rx;
w  = v*R;

% Convert units to mm^{-1}
w  = w/1e3;
kx = kx/1e3;
ky = ky/1e3;

if isfield(options,'logscale') && options.logscale
    imagesc(kx,ky,10*log10(abs(I.')/Imax))
    clim([-30 0])
else
    imagesc(kx,ky,abs(I.'))
end
hold on
plot(w(:,1),w(:,2),'r.')
xlabel('k_x (mm^{-1})')
ylabel('k_y (mm^{-1})')
drawnow

end