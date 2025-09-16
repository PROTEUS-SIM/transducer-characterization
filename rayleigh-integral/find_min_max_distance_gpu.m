function [rmin,rmax] = find_min_max_distance_gpu(source,sensor,...
    batchMemory,useGPU)
%FIND_MIN_MAX_DISTANCE_GPU finds the minimum and maximum seperation between
%source and sensor points using batched GPU computation.
%
% batchMemory is the maximum allowable memory to be allocated in bytes.
%
% useGPU: boolean, true for GPU use, false for CPU use
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

% Assign batches:
bytesPerDouble = 8;
batchSize = floor(batchMemory/(Nsource*Ndim*bytesPerDouble));
batchSize = max(1,batchSize);
i1 = 1:batchSize:Nsensor; % Start index for each batch
i2 = i1 + batchSize - 1;  % End index for each batch
i2(end) = Nsensor;
Nbatch = length(i1);      % Number of batches

% Compute pairwise distance matrix (Nsensor-by-1-by-Nsource):
r1 = reshape(source.points,[1 Nsource Ndim]);
r2 = reshape(sensor.points,[Nsensor 1 Ndim]);

if useGPU
    r1 = gpuArray(r1);
    r2 = gpuArray(r2);
end

rmin = Inf;
rmax = -Inf;

fprintf('Computing pairwise distances:     ')
for n = 1:Nbatch
    r = vecnorm(r1-r2(i1(n):i2(n),:,:),2,DIM);
    rmin = min(rmin,min(r(:)));
    rmax = max(rmax,max(r(:)));
    fprintf('\b\b\b\b%3.0f%%', n/Nbatch*100)
end
fprintf('\b\b\b\bDone\n')

rmin = gather(rmin);
rmax = gather(rmax);


end