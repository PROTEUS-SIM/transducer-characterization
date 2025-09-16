function Transducer = update_transducer(...
    Transducer,Medium,Grid,SimulationParameters)
% Divide the surfaces of the transducer elements into integration elements
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% This function assumes that the grid spacing is equal in all directions.
% Confirm that the grid is isotropic to prevent unexpected results.
if not(isequal(Grid.dx,Grid.dy,Grid.dz))
    error('Grid is not isotropic.')
end

Transducer.IntegrationDensity = SimulationParameters.IntegrationDensity;
Transducer.OnGrid = SimulationParameters.TransducerOnGrid;
Transducer = get_transducer_integration_points(Transducer, Grid);
Transducer = get_transducer_integration_delays(Transducer,Medium);

% By default, the function get_transducer_integration_points creates a
% transducer with the propagation axis along the x axis. Permute the
% coordinates, if another propagation axis is requested:
if     isfield(Transducer,'Axis') && Transducer.Axis == 2
    Transducer.integration_points = ...
        Transducer.integration_points(:,:,[3 1 2]);
elseif isfield(Transducer,'Axis') && Transducer.Axis == 3
    Transducer.integration_points = ...
        Transducer.integration_points(:,:,[2 3 1]);
end

end