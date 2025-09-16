function [S, frequencies] = define_source_grid(Transducer, Transmit, ...
    Medium, Grid, frequencies, useGPU)
%DEFINE_SOURCE_GRID converts GUI output parameter structs into a source
%that is defined as a complex distribution on a cartesian grid.
%
% S = DEFINE_SOURCE_GRID(Transducer, Transmit, Medium, Grid) returns S, a
% complex valued matrix of size Grid.Ny by Grid.Nz of which the values
% represent the complex velocity out each grid point.
%
% The source S is a band-limited representation of the transducer. This
% function currently only works for transducers defined on the plane 
% Grid.x = 0.
%
% See also define_source_frequency
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Continuous-space source representation (point cloud)
source1 = define_source_frequency(Transducer, Transmit, Medium, ...
    Grid, frequencies, useGPU);

% Grid-based source representation and weights relating points and grid
% representations
[source2, sensor_weights] = define_sensor_transducer(Transducer, Grid);
source_weights = transpose(sensor_weights);

frequencies = source1.frequencies;
source2.mask = repmat(source2.mask,[1 1 1 length(frequencies)]);

% Project the continuous-space source representation onto the grid
if strcmp(source1.type,'velocity')
    A = source1.velocities*Transducer.integration_weights;
else
    A = source1.pressures*Transducer.integration_weights;
end
S = double(source2.mask);
S(source2.mask) = gather(full(source_weights*A));
S = squeeze(S);

end