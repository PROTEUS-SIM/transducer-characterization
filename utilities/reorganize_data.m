% Reorganise the experimental data for use in this repository.
%
% This script is specific to the experimental setup that was used for
% Chapter 6 of my thesis, see Section 6.2.5 Data acquisition.
%
% This script is aimed at colleagues working at the Physics of Fluids
% group, but may be useful to other users as well, as it shows how
% experimental data should be organised in order to be used in this
% code repository.
%
% The data used in this script is not part of this repository and is
% available upon request.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear
clc

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
addpath(rootDirectory); clear currentDirectory
PATHS = path_setup_calibration();

saveNameScanFull = fullfile(PATHS.Data,'FullScanData.mat');
saveNameScan     = fullfile(PATHS.Data,'ScanData.mat');
saveNameTransmit = fullfile(PATHS.Data,'TransmitData.mat');
saveNameReceive  = fullfile(PATHS.Data,'ReceiveData.mat');

% Folder containing the oscilloscope data from the hydrophone scan:
hydrophoneDataDir = fullfile(fileparts(rootDirectory),...
    'ImpulseResponseTransducerSurface',...
    'ImpulseResponseTransducerSurface');

% Verasonics data transmit:
loadNameTransmit = fullfile(fileparts(rootDirectory),...
    'ImpulseResponseTransducerSurface',...
    'Nathan_20230325',...
    'SetUpP4_1_ImpulseResponseTriggeredVSX.mat');

% Verasonics data receive:
loadNameReceive = fullfile(fileparts(rootDirectory),...
    'ImpulseResponseExperiments',...
    'Nathan_20221220',...
    'ImpulseResponse2wy.mat');

% =========================================================================
% VERASONICS DATA
% =========================================================================

% Transmit experiment data
load(loadNameTransmit,'Resource','TPC','Trans','TW')

% Verasonics system clock speed:
Transmit.SamplingRate = Resource.VDAS.sysClk*1e6;

% Driving pulse [V*s]:
V_offset = -0.6; % Offset Verasonics transmit voltage
voltage = TPC(1).hv + V_offset;
tx = 1;
Transmit.VoltageSignal = voltage*TW(tx).TriLvlWvfm;

save(saveNameTransmit,'Trans','Transmit')

%--------------------------------------------------------------------------

% Receive experiment data
% Also includes a transmit signal, as the transmit voltage was different
% for this experiment.

load(loadNameReceive,'Resource','TPC','Trans','TW')

% Verasonics system clock speed:
Transmit.SamplingRate = Resource.VDAS.sysClk*1e6;

% Driving pulse [V*s]:
voltage = TPC(1).hv + V_offset;
tx = 1;
Transmit.VoltageSignal = voltage*TW(tx).TriLvlWvfm;

[Receive.RF,Receive.SamplingRate] = load_receive_data(loadNameReceive);

save(saveNameReceive,'Trans','Transmit','Receive')


% =========================================================================
% HYDROPHONE DATA
% =========================================================================

% Load the scan parameters:
load(fullfile(hydrophoneDataDir,'Parameters.mat'))

% Data file format: 'pos_<Xidx>_<Yidx>_<Zidx>_measurementData.mat'
% With Xidx, Yidx, Zidx the coordinate indices.

% List of all the oscilloscope data files in the hydrophone data directory:
filelist = dir(fullfile(hydrophoneDataDir, 'pos_*'));

% Size of the scan:
Nx = ScanProperties.size(1);
Ny = ScanProperties.size(2);
Nz = ScanProperties.size(3);

% Number of time samples per RF signal:
Nt = ScopeInfo.NumberOfPostTriggerSamples;

% Number of waveforms:
Nw = ScopeInfo.NumberOfTriggers/ScopeInfo.NumberOfSequenceRepetitions;

% Preallocate matrix for holding voltages:
V = zeros(Nx,Ny,Nz,Nt,Nw);

for fileIdx = 1:length(filelist)

    disp([num2str(fileIdx) '/' num2str(length(filelist)) ': Processing file ' filelist(fileIdx).name])
    load(fullfile(filelist(fileIdx).folder,filelist(fileIdx).name))

    filenameParts = split(filelist(fileIdx).name,'_');
    Xidx = str2double(filenameParts{2});
    Yidx = str2double(filenameParts{3});
    Zidx = str2double(filenameParts{4});

    V(Xidx,Yidx,Zidx,:,:) = chA;

end

sz          = ScanProperties.size;
increments  = ScanProperties.increments;
origin      = ScanProperties.origin;
minPosition = ScanProperties.minPosition;

% Scan coordinates in steps:
x = (0:(sz(1)-1))*increments(1) + minPosition(1) - origin(1);
y = (0:(sz(2)-1))*increments(2) + minPosition(2) - origin(2);
z = (0:(sz(3)-1))*increments(3) + minPosition(3) - origin(3);

% Scan coordinates in meters:
x = x*ScanProperties.MillimeterPerStep(1);
y = y*ScanProperties.MillimeterPerStep(2);
z = z*ScanProperties.MillimeterPerStep(3);

% Sampling rate
Fs = 1/ScopeInfo.timeIntervalNanoSeconds*1e9;

dx = ScanProperties.increments(1).*ScanProperties.MillimeterPerStep(1);
dy = ScanProperties.increments(2).*ScanProperties.MillimeterPerStep(2);
dz = ScanProperties.increments(3).*ScanProperties.MillimeterPerStep(3);

dx = dx*1e-3;
dy = dy*1e-3;
dz = dz*1e-3;

Grid.x = x*1e-3; Grid.Nx = size(V,1); Grid.dx = dx;
Grid.y = y*1e-3; Grid.Ny = size(V,2); Grid.dy = dy;
Grid.z = z*1e-3; Grid.Nz = size(V,3); Grid.dz = dz;

Grid = get_kspace(Grid);

Grid.Nt = size(V,4);
Grid.dt = 1/double(Fs);

medium.sound_speed = 1481; % Speed of sound in the medium [m/s];
medium.density     = 998;  % Density of the medium [kg/m^3]

% Polarity of the transmit pulses:
polarity = [1 -1];

% Save full data set:
save(saveNameScanFull,'Fs','Grid','V','medium','polarity','x','y','z')


%--------------------------------------------------------------------------
% Save streamlined data set for publication
%--------------------------------------------------------------------------

% Average positive and negative polarity pulses:
V = (V(:,:,:,:,1)*polarity(1) + V(:,:,:,:,2)*polarity(2))/2;

% Remove dimensions of length 1
voltage_mV = squeeze(V);

% Rotate coordinate system such that z aligns with the propagation axis:
Grid.full_size = [Grid.Nx; Grid.Ny; Grid.Nz];
Grid = permute_grid(Grid,[2 3 1]);
Grid = get_kspace(Grid);

% Save data:
save(saveNameScan,'Grid','medium','voltage_mV')


% Notes for personal archive Nathan Blanken:
%
% Original location of scan data output:
% ImpulseResponseTransducerSurface/FullScanData.mat
%
% Original location of the filter file:
% ImpulseResponseTransducerSurface/BLP-15+.txt
%
% Original location of data/PROTEUS-I/PA_XZ.mat
% PROTEUS/virtual-transducer/P4-1_MATLAB-files-2/PA_XZ.mat
%
% Original location of data/PROTEUS-I/pressure_maps.mat
% PROTEUS/virtual-transducer/PressureField_2cycles/pressure_maps.mat
