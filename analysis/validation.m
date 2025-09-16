% Validation of the virtual transducer on independent experimental data
% (Section 6.3.3 of my thesis).
%
% Creates the base plots for Fig. 6.8 of my thesis.
%
% Dependencies:
% data/PROTEUS-I
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%#ok<*UNRCH>

clear
clc
close all

saveFigure = false;
rerun = true;

% Simulation B was performed with a lower spatial resolution than
% simulation A, resulting in a slightly lower output pressure. The
% undersampling factor compensates for this.
undersamplingfactor = 1.054;
fitfactor = 1.14;
factorA = 1;
factorB = undersamplingfactor*fitfactor;

waveformIndex = 17;
centreLineIndexExperiment = 17;
        
% Maximum pressure to display in kPa:
maxPressure_kPa = 500;

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
addpath(rootDirectory); clear rootDirectory
PATHS = path_setup_calibration();

savefolder = 'pressure-maps-PROTEUS';

if rerun
    % Run PROTEUS pressure field simulation with the new calibration
    % results
    settingsfile = fullfile(currentDirectory,'GUI_output_parameters.mat');
    main_pressure_field(settingsfile, savefolder);

    % Remove any existing data
    if exist(savefolder,'dir')
        delete('pressure-maps-PROTEUS/*')
        rmdir('pressure-maps-PROTEUS')
    end

    % Move the results from the PROTEUS results folder to the current
    % directory
    addpath(PATHS.PROTEUS)
    PATHS_PROTEUS = path_setup();
    movefile(fullfile(PATHS_PROTEUS.ResultsPath,savefolder),...
        currentDirectory)
    clear currentDirectory PATHS_PROTEUS
    rmpath(PATHS.PROTEUS)
end

ExperimentFile  = fullfile(PATHS.Data,'PROTEUS-I','PA_XZ.mat');
SimulationFileA = fullfile(PATHS.Data,'PROTEUS-I','pressure_maps.mat');
SimulationFileB = fullfile(savefolder,'pressure_maps.mat');

% Load the experimental and simulation data:
load(ExperimentFile)
S = load(SimulationFileA);
GridA = S.Grid;
sensor_data_A_xy = S.sensor_data_xy;
sensor_data_A_xz = S.sensor_data_xz;
clear S

S = load(SimulationFileB);
GridB = S.Grid;
sensor_data_B_xy = S.sensor_data_xy;
sensor_data_B_xz = S.sensor_data_xz;
clear S

centreLineIndexSimulationA = find(GridA.z==0);
centreLineIndexSimulationB = find(GridB.z==0);

% Get peak to peak pressure map in the requested plane in kPa:
pressureMapExp = squeeze(PA_pk2pk(:,1,:,waveformIndex));
pressureMapSimA = sensor_data_A_xz.p_max - sensor_data_A_xz.p_min;
pressureMapSimB = sensor_data_B_xz.p_max - sensor_data_B_xz.p_min;

pressureMapSimA = pressureMapSimA*factorA;
pressureMapSimB = pressureMapSimB*factorB;

pressureMapSimA = pressureMapSimA/1e3; % Convert to kPa
pressureMapSimB = pressureMapSimB/1e3; % Convert to kPa

% Pressure map cross sections:
pressureProfileExp  = pressureMapExp(:, centreLineIndexExperiment);
pressureProfileSimA = pressureMapSimA(:,centreLineIndexSimulationA);
pressureProfileSimB = pressureMapSimB(:,centreLineIndexSimulationB);

% Rotate pressure maps
pressureMapSimA = transpose(pressureMapSimA);
pressureMapSimB = transpose(pressureMapSimB);
pressureMapExp  = transpose(pressureMapExp);

% Get coordinates of the pressure maps in mm:
coordinateExp1 = X;
coordinateSimA1 = GridA.x*1e3;
coordinateSimB1 = GridB.x*1e3;
label1 = 'X (mm)';

coordinateExp2 = Z;
coordinateSimA2 = GridA.z*1e3;
coordinateSimB2 = GridB.z*1e3;
label2 = 'Z (mm)';

% Subplot positions:
positionExp  = 1;
positionSimA = 2;
positionSimB = 3;

fig1 = figure(1);
subplot(3,1,positionExp)
imagesc(coordinateExp1, coordinateExp2, pressureMapExp)
clim([0 maxPressure_kPa])
axis equal
ax1 = gca;

xlim([min(coordinateExp1) max(coordinateExp1)])
ylim([min(coordinateExp2) max(coordinateExp2)])
xlabel(label1)
ylabel(label2)

subplot(3,1,positionSimA)
imagesc(coordinateSimA1, coordinateSimA2, pressureMapSimA);
clim([0 maxPressure_kPa])
axis equal
xlim([min(coordinateExp1) max(coordinateExp1)])
ylim([min(coordinateExp2) max(coordinateExp2)])
xlabel(label1)
ylabel(label2)
ax2 = gca;

subplot(3,1,positionSimB)
imagesc(coordinateSimB1, coordinateSimB2, pressureMapSimB);
clim([0 maxPressure_kPa])
axis equal
xlim([min(coordinateExp1) max(coordinateExp1)])
ylim([min(coordinateExp2) max(coordinateExp2)])
xlabel(label1)
ylabel(label2)
ax3 = gca;

h = colorbar;
h.Location = 'south';
h.Label.String = 'Peak-to-peak pressure (kPa)';

figure
plot(coordinateSimA1, pressureProfileSimA)
hold on
plot(coordinateSimB1, pressureProfileSimB)
plot(coordinateExp1, pressureProfileExp)
xlabel(label1)
ylabel('Pressure (kPa)')
fig4 = gcf;
ax4 = gca;

if not(saveFigure)
    return
end

saveFolder = string(fullfile(PATHS.Figures,'fig'));

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

savefig(fig1, saveFolder + filesep + "Fig8a.fig")
savefig(fig4, saveFolder + filesep + "Fig8d.fig")
