function beta = compute_obliquity_factor(...
    source_points, sensor_points, normal, baffle)
%COMPUTE_OBLIQUITY_FACTOR computes the obliquity factor for pairs of source
%points on the transducer surface and sensor points.
%
% beta = COMPUTE_OBLIQUITY_FACTOR(source_points, sensor_points, normal,
% baffle) returns the obliquity factors beta (Nsensor-by-Nsource) for pairs
% of source points (Nsource-by-3) and sensor points (Nsensor-by-3) with
% normal (Nsource-by-3) the normal to the transducer surface at each source
% point and baffle the baffle type: 'rigid', 'soft', or 'free'.
%
% Obliquity factors derived from J.L. San Emeterio and L. G. Ullate:
% Diffraction impulse response of rectangular transducers J. Acoust. Soc.
% Am. 92 (2), Pt. 1, August 1992. Obtained by applying Eq. 5 to a
% harmonically oscillating velocity. See https://doi.org/10.1121/1.403990.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

Nsource = size(source_points,1);      % Number of source points
Nsensor = size(sensor_points,1);      % Number of sensor points
Ndim    = size(source_points,2);      % Number of spatial dimensions (3)

% Matrix dimension corresponding to the spatial coordinates:
DIM = 3;

% Cosine of the angles between lines from the integration points and the
% observation point and the normal to the transducer:
if strcmp(baffle,'soft') || strcmp(baffle,'free')
    
    % Compute pairwise distance matrix (Npoints2-by-Npoints1):
    r1 = reshape(source_points,[1,Nsource,Ndim]);
    r2 = reshape(sensor_points,[Nsensor,1,Ndim]);
    r  = vecnorm(r1-r2,2,DIM);
    
    % Normal vector pointing into the medium for each point on the
    % transducer:
    n = reshape(normal,[1,Nsource,Ndim]);
    n = repmat(n,[Nsensor,1,1]);
    
    costheta = dot(r2-r1,n,DIM)./r;
end

switch baffle
    case 'rigid'
        beta = ones(Nsensor,Nsource);
    case 'soft'
        beta = costheta;
    case 'free'
        beta = (costheta + ones(Nsensor,Nsource))/2;
end

end