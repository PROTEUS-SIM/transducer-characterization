function [thetaX,thetaY,thetaZ] = find_angles(A,thetaZ,kx,ky,kxi,kyi,k)
%FIND_ANGLES finds the rotation angles that define the rotation of the
%sensor plane with respect to the source plane.
%
% If A is the angular spectrum of sensor data obtained on a rotated sensor
% plane, [thetaX, thetaY, thetaZ] = FIND_ANGLES(A,thetaZ,kx,ky,kxi,kyi,k)
% finds the rotation angles thetaX, thetaY, and thetaZ that define the
% rotation of the sensor plane with respect to the source plane.
%
% The z axis is the propagation axis, and the x and y axis define the
% source plane.
%
% The sensor plane is obtained by a rotation of thetaX radians about the x
% axis, followed by a rotation of thetaY radians about the y axis, followed
% by a rotation of thetaZ about the z axis.
%
% kx and ky are the grid vectors of the k-space grid.
%
% k is the wavenumber: k = 2*pi/lambda, where lambda is the wavelength of
% the wave.
%
% kxi and kyi define a set of integration points (kxi,kyi). These are
% points in k-space where the intensity of the angular spectrum would be
% high if the sensor plane were not rotated. The algorithm assumes that
% (kxi,kyi)=(0,0) corresponds to the highest intensity in the angular
% spectrum.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Ensure kxi and kyi are column vectors:
if isrow(kxi)
    kxi = transpose(kxi);
end

if isrow(kyi)
    kyi = transpose(kyi);
end

% Convert the integration points to vectors in 3D k-space:
ki = [kxi kyi sqrt(k^2 -kxi.^2 - kyi.^2)*sign(k)];

%==========================================================================
% COMPUTE THETAX AND THETAY
%==========================================================================
% Find thetaX and thetaY from the highest intensity point in the angular
% spectrum.

Nx = length(kx); dkx = mean(diff(kx));
Ny = length(ky); dky = mean(diff(ky));

% Only search within the circle kx^2 + ky^2 = k^2:
[Kx,Ky] = ndgrid(kx,ky);
dk = max(dkx,dky); % Margin for filter
A(Kx.^2 + Ky.^2 >= (abs(k)-dk)^2) = 0;

[~,imax] = max(A(:));
[imax,jmax] = ind2sub([Nx,Ny],imax);

kx_max = kx(imax);
ky_max = ky(jmax);

% With R = Rz*Ry*Rx, these relations hold for a wave moving along the
% z-axis:
%
% kx_max = -k*sin(thetaY)
% ky_max =  k*sin(thetaX)*cos(thetaY)
%
% With -pi/2<= thetaY <= pi/2, these can be rewritten as:
thetaX = asin(ky_max/sqrt(k^2-kx_max^2)*sign(k));
thetaY = asin(-kx_max/k);

%==========================================================================
% COMPUTE THETAZ
%==========================================================================
% For each element of thetaZ, compute the rotation matrix, rotate the
% integration vectors, and integrate the angular spectrum over the
% integration vectors. Return the element of thetaZ for which the
% integration result is maximum.

% Get rotation matrices:
Rx = get_rotation_matrix(thetaX,1);
Ry = get_rotation_matrix(thetaY,2);

[Kx,Ky] = ndgrid(kx,ky);
F = griddedInterpolant(Kx,Ky,A);

% Array for holding integration results:
S = zeros(size(thetaZ));

for n = 1:length(thetaZ) 
    
    Rz = get_rotation_matrix(thetaZ(n),3);
    R = Rz*Ry*Rx;      % Trial rotation matrix
    ki_rotated = ki*R; % Rotated integration vectors
       
    S(n) = sum(F(ki_rotated(:,1),ki_rotated(:,2)));
    
end

[~,imax] = max(S);
thetaZ = thetaZ(imax);

end
