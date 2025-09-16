function [rmin,rmax] = find_min_max_distance(source,sensor)
%FIND_MIN_MAX_DISTANCE finds the minimum and maximum seperation between
%source and sensor points.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

Nsource = size(source.points,1);      % Number of source points
Nsensor = size(sensor.points,1);      % Number of sensor points
Ndim    = size(source.points,2);      % Number of spatial dimensions (3)

% Matrix dimension corresponding to the spatial coordinates:
DIM = 3;

% Compute pairwise distance matrix (Nsensor-by-1-by-Nsource):
r1 = reshape(source.points,[1 Nsource Ndim]);
r2 = reshape(sensor.points,[Nsensor 1 Ndim]);
r = vecnorm(r1-r2,2,DIM);
rmin = min(r(:));
rmax = max(r(:));

end