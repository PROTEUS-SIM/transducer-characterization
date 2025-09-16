% This script shows a comparison of propagation methods and source
% presentations as described in Section 6.2.3 of my thesis.
%
% The data for this comparison can be generated with
% comparison_propagation_methods.m
%
% Creates the base plots for Fig. 6.2 of my thesis.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear; clc; close all

saveFigure = false;

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
resultsDirectory = fullfile(currentDirectory,'data');
addpath(rootDirectory); clear currentDirectory rootDirectory
PATHS = path_setup_calibration;

filename1 = fullfile(resultsDirectory,'sensor_data_dipole_rayleigh.mat');
filename2 = fullfile(resultsDirectory,'sensor_data_dipole_angular.mat');
filename3 = fullfile(resultsDirectory,'sensor_data_dipole_k_Wave.mat');
filename4 = fullfile(resultsDirectory,'sensor_data_monopole_rayleigh.mat');
filename5 = fullfile(resultsDirectory,'sensor_data_monopole_angular.mat');
filename6 = fullfile(resultsDirectory,'sensor_data_monopole_k_Wave.mat');

figure
hold on

load(filename4,'t')
N = length(t);
p_monopole = zeros(3,N);
p_dipole   = zeros(3,N);

displayUnitTime = 1e-6;
displayUnitPressure = 1e3;
load(filename1,'t','p')
plot(t/displayUnitTime,p/displayUnitPressure)
p_dipole(1,:) = p;
load(filename2,'t','p')
plot(t/displayUnitTime,p/displayUnitPressure)
p_dipole(2,:) = p;
load(filename3,'t','p')
plot(t/displayUnitTime,p/displayUnitPressure)
p_dipole(3,:) = p;

load(filename4,'t','p')
plot(t/displayUnitTime,p/displayUnitPressure)
p_monopole(1,:) = p;
load(filename5,'t','p')
plot(t/displayUnitTime,p/displayUnitPressure)
p_monopole(2,:) = p;
load(filename6,'t','p')
plot(t/displayUnitTime,p/displayUnitPressure)
p_monopole(3,:) = p;

xlim([0 15])
ylim([-250 250])
xlabel('Time (microseconds)')
ylabel('Pressure (kPa)')

legend('dipole Rayleigh','dipole angular','dipole k-Wave', ...
    'monopole Rayleigh','monopole angular','monopole k-Wave')

p_monopole_diff = reshape(p_monopole,1,3,N) - reshape(p_monopole,3,1,N);
p_dipole_diff   = reshape(p_dipole,1,3,N)   - reshape(p_dipole,3,1,N);

disp('Maximum difference')
disp(max(p_monopole_diff(:))/max(p_monopole(:)))
disp(max(p_dipole_diff(:))/max(p_dipole(:)))

if not(saveFigure)
    return
end

saveFolder = string(fullfile(PATHS.Figures,'fig'));

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

savefig(gcf, saveFolder + filesep + "Fig2.fig")