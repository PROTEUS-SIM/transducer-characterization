function PATHS = path_setup_calibration(varargin)
% Load paths to modules and data files
%
% For installation: modify the variable PATHS.AcousticModule to correspond
% to the full path on your system.
%
% The paths PATHS.kWave and PATHS.PROTEUS are optional and can be left
% empty ('').
%
% Running PATHS = path_setup_calibration() or PATHS =
% path_setup_calibration('addpath') adds the required folders to the MATLAB
% path.
%
% Running PATHS = path_setup_calibration('rmpath') removes these folders
% from the MATLAB path again.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Modules (required), added to MATLAB path
PATHS.AcousticModule   = '/home/user/PROTEUS/acoustic-module';

% Modules (optional), not added to MATLAB path
PATHS.kWave            = '/home/user/k-wave-toolbox-version-1.3/k-Wave';
PATHS.PROTEUS          = '/home/user/PROTEUS';

% =========================================================================
% ---------------------------No need to modify-----------------------------
% =========================================================================

% Full path to current function:
currentFile = mfilename('fullpath');

% Full path to the start directory:
rootDirectory = fileparts(currentFile);

% Set up data path:
PATHS.Data = fullfile(rootDirectory,'data');
if not(isfolder(PATHS.Data))
    mkdir(PATHS.Data)
end

% Set up results path:
PATHS.Results = fullfile(rootDirectory,'results');
if not(isfolder(PATHS.Results))
    mkdir(PATHS.Results)
end

% Set up checkpoints path:
PATHS.Checkpoints = fullfile(rootDirectory,'checkpoints');
if not(isfolder(PATHS.Checkpoints))
    mkdir(PATHS.Checkpoints)
end

% Figure path:
PATHS.Figures = fullfile(rootDirectory,'figures');

% Add folders to or remove folders from MATLAB path
if isempty(varargin) || strcmp(varargin,'addpath')
    addpath(...
        rootDirectory, ...
        fullfile(rootDirectory,'angular-spectrum'), ...
        fullfile(rootDirectory,'rayleigh-integral'), ...
        fullfile(rootDirectory,'utilities'), ...
        PATHS.AcousticModule)
elseif strcmp(varargin,'rmpath')
    rmpath(...
        rootDirectory, ...
        fullfile(rootDirectory,'angular-spectrum'), ...
        fullfile(rootDirectory,'rayleigh-integral'), ...
        fullfile(rootDirectory,'utilities'), ...
        PATHS.AcousticModule)
else
    error('Invalid input argument.')
end

end

