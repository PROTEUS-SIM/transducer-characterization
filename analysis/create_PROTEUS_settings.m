% Write a PROTEUS settings file based on the transducer characterization
% results that can be used in the script validation.m
%
% The settings file can also be created entirely via the graphical user
% interface MainGUI in the PROTEUS repository.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear; clc; close all

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
addpath(rootDirectory); clear rootDirectory
PATHS = path_setup_calibration();

addpath(PATHS.PROTEUS)
PATHS_PROTEUS = path_setup();
addpath(PATHS_PROTEUS.GUIfunctions)

% -------------------------------------------------------------------------
% Acquisition and Microbubble
% -------------------------------------------------------------------------

% Default properties for Acquisition and Microbubble, not necessary for
% pressure field simulation
Acquisition   = reset_acquisition();
Microbubble   = reset_microbubble();

% -------------------------------------------------------------------------
% Medium
% -------------------------------------------------------------------------

% Select water
Medium = reset_medium();
Medium.Tissue = 'Water';
Medium = assign_medium_properties(Medium);

% Make the medium fully homogeneous
Medium.Inhomogeneity = 0;
Vessel.Tissue = 'Water';
Vessel = assign_medium_properties(Vessel);
Vessel = assign_liquid_properties(Vessel,Vessel.Tissue);
Medium.Vessel = Vessel; clear Vessel
Medium.SpeedOfSoundMinimum = Medium.SpeedOfSound;
Medium.SpeedOfSoundMaximum = Medium.SpeedOfSound;

Medium.Save = false;

% -------------------------------------------------------------------------
% Transducer
% -------------------------------------------------------------------------

T = load(fullfile(PATHS.Results,'Transducer.mat'),'Transducer');
Transducer                  = reset_transducer();
Transducer.Type             = 'Custom linear';
Transducer.NumberOfElements = T.Transducer.NumberOfElements;
Transducer.Pitch            = T.Transducer.Pitch;
Transducer.ElementWidth     = T.Transducer.ElementWidth;
Transducer.ElementHeight    = T.Transducer.ElementHeight;
Transducer.ElevationFocus   = T.Transducer.ElevationFocus;

Transducer.ImpulseResponseType = 'Load file';
load(fullfile(PATHS.Results,'IR_receive_processed.mat'),'IR')
Transducer.ReceiveImpulseResponse = IR;
load(fullfile(PATHS.Results,'IR_transmit_processed.mat'),'IR')
Transducer.TransmitImpulseResponse = IR;

% -------------------------------------------------------------------------
% Transmit
% -------------------------------------------------------------------------

% Transmit data used in PROTEUS Part I:
transmitdatafolder = fullfile(PATHS.PROTEUS,'example_custom_data');
apodizationFile = fullfile(transmitdatafolder,'custom_apodization.mat');
delaysFile      = fullfile(transmitdatafolder,'custom_delays.mat');

Transmit = reset_transmit(Transducer);
Transmit.Advanced = true;
Transmit.AmplitudeMode = 'Voltage';
Transmit.VoltageAmplitude = 6.9;
Transmit.ApodizationType = 'Load apodization';
Transmit.DelayType       = 'Load delays';
load(apodizationFile,'apodization');
load(delaysFile,'delays')
Transmit.Apodization = apodization; clear apodization
Transmit.Delays = delays; clear delays
Transmit = get_voltage_signal(Transmit);
Transmit = get_pressure_signal(Transmit,Transducer);

% -------------------------------------------------------------------------
% Geometry
% -------------------------------------------------------------------------

Geometry = reset_geometry([],PATHS_PROTEUS);
Geometry.startDepth = 0.115;
Geometry.Domain.Margin = 2*Medium.SpeedOfSound/Transmit.CenterFrequency;
Geometry = compute_simulation_domain(Geometry, Transducer, Transmit);

% -------------------------------------------------------------------------
% Simulation Parameters
% -------------------------------------------------------------------------

SimulationParameters = reset_simulation_parameters();
SimulationParameters.PointsPerWavelength = 4;
SimulationParameters.Solver = '3DG';
SimulationParameters = update_simulation_parameters(...
    SimulationParameters, Medium, Transmit.CenterFrequency);

% -------------------------------------------------------------------------
% Write file
% -------------------------------------------------------------------------

filename = 'GUI_output_parameters.mat';
write_to_file(Microbubble,SimulationParameters,Geometry,...
                Transducer,Acquisition,Medium,Transmit,[],filename)

rmpath(PATHS.PROTEUS)
rmpath(PATHS_PROTEUS.GUIfunctions)
