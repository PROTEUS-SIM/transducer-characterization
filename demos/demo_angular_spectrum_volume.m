% This script demonstrates how the angular spectrum method can be used to
% rapidly generate the volume of the transducer output field by propagating
% to closely spaced planes.
%
% This result is not included in my thesis.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear; clc

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
addpath(rootDirectory); clear currentDirectory rootDirectory
PATHS = path_setup_characterization();

load('GUI_output_parameters.mat', 'Geometry', ...
    'Medium', 'SimulationParameters', 'Transducer', 'Transmit')

expansion_factor = 2; % Grid expansion factor
Geometry.Axis = 3; % Propagation axis
Grid = define_grid_demos(SimulationParameters,Geometry,expansion_factor);

Transducer.Axis = 3; % Propagation axis
Transducer=update_transducer(Transducer,Medium,Grid,SimulationParameters);

colorMap = hot;
colorMap = flipud(colorMap);

%==========================================================================
% SOURCE REPRESENTATION
%==========================================================================

f0 = Transmit.CenterFrequency;
Transmit.ContinuousWave = true;
S = define_source_grid(Transducer, Transmit, Medium, Grid, f0, false);

%==========================================================================
% ANGULAR SPECTRUM METHOD
%==========================================================================

c0     = Medium.SpeedOfSound;
rho    = Medium.Density;

% Compute the wavenumber of propagation:
k  = 2*pi*f0/c0;

Grid = get_kspace(Grid);
Nx = Grid.Nx; dkx = Grid.dkx;
Ny = Grid.Ny; dky = Grid.dky;

%==========================================================================
% SHOW RESULTS
%==========================================================================

displayUnit = 1e-3; % [m]

d = 0:0.005:0.15;

G1 = compute_spectral_propagator(Nx,dkx,Ny,dky,k,d,true,'vp');

% Apply forward propagation:
S1 = angular_spectrum(S,Grid.x,Grid.y);
S1 = S1.*G1*rho*c0;
S1 = angular_spectrum_inverse(S1,Grid.x,Grid.y);

figure

for n = 1:length(d)
    
    % Plot the pressure field at the transducer surface:
    [X,Y,Z] = ndgrid(Grid.x,Grid.y,d(n));
    X = X/displayUnit;
    Y = Y/displayUnit;
    Z = Z/displayUnit;
    
    alphaChannel = abs(S1(:,:,n));
    alphaChannel = alphaChannel/max(alphaChannel(:));
    p = surf(Z,X,Y,abs(S1(:,:,n)),...
        'AlphaData',alphaChannel,'FaceAlpha','flat');
    p.EdgeColor = 'none';

    hold on
    
end

hold off
colormap(colorMap)

xlabel('z (mm)')
ylabel('x (mm)')
zlabel('y (mm)')

ax = gca;
set(ax,'DataAspectRatio',[1 1 1])
