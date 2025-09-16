% Demonstration of the angular spectrum method.
%
% Creates the base plots for Fig. 6.1 of my thesis.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%#ok<*UNRCH>

clear; clc; close all

saveFigure = false;

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
addpath(rootDirectory); clear currentDirectory rootDirectory
PATHS = path_setup_characterization();

if saveFigure == true
    d1 = 20e-3;  % Propagation distance 1 [m]
    d2 = -20e-3; % Propagation distance 2 [m]
    
    ppwl = 6; % New number of points per wavelength
    expansion_factor = 4; % Grid expansion factor
else
    d1 = 3e-3;   % Propagation distance 1 [m]
    d2 = -3e-3;  % Propagation distance 2 [m]
    
    ppwl = 6; % New number of points per wavelength
    expansion_factor = 2; % Grid expansion factor
end

load('GUI_output_parameters.mat', 'Geometry', ...
    'Medium', 'SimulationParameters', 'Transducer', 'Transmit')

Geometry.Axis = 3; % Propagation axis
Grid = define_grid_demos(SimulationParameters,Geometry,...
    expansion_factor,ppwl);

Transducer.Axis = 3; % Propagation axis
Transducer=update_transducer(Transducer,Medium,Grid,SimulationParameters);

%==========================================================================
% SOURCE REPRESENTATION
%==========================================================================

f0 = Transmit.CenterFrequency;
Transmit.ContinuousWave = true;
S = define_source_grid(Transducer, Transmit, Medium, Grid, f0, false);

%==========================================================================
% ANGULAR SPECTRUM METHOD
%==========================================================================

f0     = Transmit.CenterFrequency;
c0     = Medium.SpeedOfSound;
rho    = Medium.Density;

% Compute the wavenumber of propagation:
k  = 2*pi*f0/c0;

Grid = get_kspace(Grid);
Nx = Grid.Nx; kx = Grid.kx; dkx = Grid.dkx;
Ny = Grid.Ny; ky = Grid.ky; dky = Grid.dky;

% Forward propagator:
G1 = compute_spectral_propagator(Nx,dkx,Ny,dky,k,d1,true,'vp');

% Backward propagator:
G2 = compute_spectral_propagator(Nx,dkx,Ny,dky,k,d2,true,'pv');

% Filter evanescent waves for stability:
G2 = filter_evanescent_waves(G2,Grid,k);

% Apply forward propagation:
Sfft = angular_spectrum(S,Grid.x,Grid.y);
Sfft2 = Sfft.*G1*rho*c0;
S2 = angular_spectrum_inverse(Sfft2,Grid.x,Grid.y);

% Apply backward propagation:
Sfft3 = angular_spectrum(S2,Grid.x,Grid.y);
Sfft3 = Sfft3.*G2/(rho*c0);
S3 = angular_spectrum_inverse(Sfft3,Grid.x,Grid.y);


%==========================================================================
% SHOW RESULTS
%==========================================================================

dynamic_range_dB = 45;
Smax = max(abs(S(:)));

fig1 = figure;
imagesc(Grid.x*1e3,Grid.y*1e3,abs(S)');
clim([0 Smax])
xlabel('Lateral, x (mm)')
ylabel('Elevation, y (mm)')
title('Original source')
axis equal
axis tight

fig4 = figure;
S0 = max(abs(Sfft(:)));
imagesc(kx*1e-3, ky*1e-3, 20*log10(abs(Sfft)/S0)');
xlabel('k_x (mm^{-1})')
ylabel('k_y (mm^{-1})')
title('Spectrum of source')
colorbar
clim([-dynamic_range_dB 0])
hold on
theta = (0:1:360)*2*pi/360;
plot(k*cos(theta)*1e-3,k*sin(theta)*1e-3,'w');
axis equal
axis tight

fig5 = figure;
S0 = max(abs(Sfft2(:)));
imagesc(kx*1e-3, ky*1e-3, 20*log10(abs(Sfft2)/S0)');
xlabel('k_x (mm^{-1})')
ylabel('k_y (mm^{-1})')
title('Spectrum after applying forward spectral propagator')
colorbar
clim([-dynamic_range_dB 0])
axis equal
axis tight

fig2 = figure;
imagesc(Grid.x*1e3,Grid.y*1e3,abs(S2)');
xlabel('Lateral, x (mm)')
ylabel('Elevation, y (mm)')
title('Propagated field')
axis equal
axis tight

fig6 = figure;
S0 = max(abs(Sfft3(:)));
imagesc(kx*1e-3, ky*1e-3, 20*log10(abs(Sfft3)/S0)');
xlabel('k_x (mm^{-1})')
ylabel('k_y (mm^{-1})')
title('Spectrum after applying backward spectral propagator')
colorbar
clim([-dynamic_range_dB 0])
axis equal
axis tight

fig3 = figure;
imagesc(Grid.x*1e3,Grid.y*1e3,abs(S3)');
clim([0 Smax])
xlabel('Lateral, x (mm)')
ylabel('Elevation, y (mm)')
title('Reconstructed source through backpropagation')
axis equal
axis tight

if not(saveFigure)
    return
end

saveFolder = string(fullfile(PATHS.Figures,'fig'));

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

savefig(fig1, saveFolder + filesep + "Fig1a.fig")
savefig(fig2, saveFolder + filesep + "Fig1b.fig")
savefig(fig3, saveFolder + filesep + "Fig1c.fig")
savefig(fig4, saveFolder + filesep + "Fig1d.fig")
savefig(fig5, saveFolder + filesep + "Fig1e.fig")
savefig(fig6, saveFolder + filesep + "Fig1f.fig")