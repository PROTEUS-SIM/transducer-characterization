function P = compute_average_receive_pressure(...
    P,Grid,medium,Transducer,options)
% This function computes the simulated receive pressure after reflection by
% the metal plate. This corresponds to the following steps in the blue part
% of Fig. 4 of the arXiv preprint:
% - Reverse lens delays
% - Compute average receive pressure
%
% The output data corresponds to Fig. 7d of the arXiv preprint.
%
% input
% - P: time-domain pressure data in sensor plane (Nx-by-Ny-by-Nsamples)
%
% output
% - P: average pressure per transducer element (Nchannels-by-Nsamples)
%
% This function leverages tools from PROTEUS to set up a Transducer
% structure. This structure can be used to map pressure data on a uniform
% grid onto the transducer elements.
%
% NOTE: This script depends on an updated version of compute_RF, which
% exists in the branch transducer-calibration of PROTEUS.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

N = size(P,3);

run_param.RF_type = 'pressure';
if options.useGPU
    run_param.DATA_CAST_RF = 'gpuArray-double';
else
    run_param.DATA_CAST_RF = 'double';
end

% Simulation settings in PROTEUS format
Medium.SpeedOfSound = medium.sound_speed;
Medium.Density      = medium.density;
SimulationParameters.IntegrationDensity = 1.5;
SimulationParameters.TransducerOnGrid = false;
Grid.dz = Grid.dx; % Make grid spacing identical in all directions

% Set up a fully defined PROTEUS transducer
Transducer.Axis = 3; % Propagation axis
Transducer = update_transducer(...
    Transducer,Medium,Grid,SimulationParameters);

% Compute the mapping between sensor points on the regular grid and sensor
% points in the transducer structure:
[sensor, sensor_weights] = define_sensor_transducer(Transducer, Grid);
mask = squeeze(sensor.mask);
mask = repmat(mask,1,1,N);

% Map sensor data onto transducer elements, reverse lens delays, and apply
% spatial averaging:
sensor_data.p = reshape(P(mask),[],N);
P = compute_RF_data(Transducer,sensor_data,sensor_weights,Grid,run_param);

end
