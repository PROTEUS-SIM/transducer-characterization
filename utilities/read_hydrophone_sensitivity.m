%READ_HYDROPHONE_SENSITIVITY reads the sensitivity data from a Precision
%Acoustics fibre-optic hydrophone sensitivity file and saves the data as a
%MAT file
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear; clc

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
addpath(rootDirectory); clear currentDirectory
PATHS = path_setup_characterization();

% Location of the Precision Acoustics sensitivity data:
filename = fullfile(fileparts(rootDirectory),...
    'ImpulseResponseTransducerSurface',...
    'FOH54 + FP196-08T Sensitivity 13June19.txt');

savename = fullfile(PATHS.Data,'sensitivity.mat');

% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 3);

% Specify range and delimiter
opts.DataLines = [82, Inf];
opts.Delimiter = "\t";

% Specify column names and types
opts.VariableNames = ["Frequency", "Sensitivity", "dS"];
opts.VariableTypes = ["double", "double", "double"];

% Import the data
tbl = readtable(filename, opts);

% Convert to output type
Frequency   = tbl.Frequency;
Sensitivity = tbl.Sensitivity;

save(savename,'Frequency','Sensitivity')
