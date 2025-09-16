function memoryEstimate = get_memory_estimate(source,sensor)
%GET_MEMORY_ESTIMATE computes the memory required to compute the pairwise
%distance matrix for sensor and source points and the memory required to
%compute the integration kernel. GET_MEMORY_ESTIMATE returns the higher
%value of the two. This value is a lower bound on the required memory for
%computing the Rayleigh integral.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

Nsource = size(source.points,1);      % Number of source points
Nsensor = size(sensor.points,1);      % Number of sensor points
Ndim    = size(source.points,2);      % Number of spatial dimensions (3)
Nfreq   = size(source.frequencies,2); % Number of frequencies

bytesPerDouble = 8;
bytesPerComplexDouble = 16;

% Memory estimate for the pairwise distance computation:
memoryEstimate1 = Nsensor*Nsource*Ndim*bytesPerDouble;

% Memory estimate for the kernel computation:
memoryEstimate2 = Nsensor*Nsource*Nfreq*bytesPerComplexDouble;

memoryEstimate = max(memoryEstimate1,memoryEstimate2);

end