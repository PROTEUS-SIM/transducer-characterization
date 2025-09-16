%MAIN_PIPELINE - A full transducer characterization pipeline from
%hydrophone data to virtual transducer
%
% This script is a demonstration of the full characterization pipeline
% shown in Fig. 6.4 of my thesis. Running the script reproduces the results
% shown in Section 6.3.
%
% It is possible (and recommended) to run only a part of this pipeline,
% save the data, and resume later through the use of checkpoints. See the
% relevant section in the script for more information on how to use them.
%
% All data are in SI units, unless specified otherwise.
%
% Radio-frequency data:
% - P: pressure data derived from hydrophone data
% - V: transducer surface normal velocity derived from hydrophone data
% - y: electronic transmit/receive RF data from ultrasound scanner with
%      sampling rate Fs
% - IR: transmit or receive impulse response
%
% Model parameters:
% - Transducer: struct holding transducer parameters in PROTEUS format
% See: https://github.com/PROTEUS-SIM/PROTEUS/blob/main/documentation/
% TransducerGUI.md
%
% Measurement parameters:
% - Grid: struct holding the grid properties of the hydrophone scan plane,
%   in PROTEUS format. Grid.dt corresponds to the sampling interval of the
%   hydrophone data. See also: define_grid
% - ScanPlane: struct holding the orientation of the hydrophone scan plane,
%   with fields thetaX, thetaY, thetaZ, x0, y0, z0
% - VirtualReceiver: struct holding the orientation of the virtual receiver
%   plane, with fields thetaX, thetaY, thetaZ, x0, y0, z0
% - medium: struct holding the properties of the propagation medium, with
%   fields density and sound_speed
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear
clc
close all

% Number of iterations of backward propagation and transducer model
% parameter estimation:
N_iterations = 2;

PATHS = path_setup_characterization();

showFigures = true;
saveFigures = false;

options.useGPU = true;
options.showFigure = showFigures;

% Frequency component to show in plots of three-dimensional data
options.plotFrequency = 2.5e6;

%==========================================================================
% Checkpoints
%==========================================================================
% Checkpoints can be used to run only a part of the pipeline.
%
% The meaning of the checkpoint values:
% 0: terminate pipeline at checkpoint
% 1: run code at checkpoint
% 2: load data from previous run at checkpoint
%
% The number of elements in checkpoint2 and checkpoint3 should equal
% N_iterations for each.
%
% Note that checkpoint2(1) is followed by checkpoint3(1), checkpoint2(2),
% and checkpoint3(2) respectively.
%
% It is recommended to start with checkpoint1 set to value 1, and all
% others to value 0 and then raise checkpoint1 to value 2 and checkpoint2
% to value 1, etc. This will provide the most interactive experience.
%
% It is not recommended to run the full pipeline (all checkpoints set to 1)
% while saveFigures = true, as this will generate a massive number of
% figures at once.
%
% Completion time:
%
% Checkpoint 1 takes on the order of a minute.
%
% Checkpoints 2 and checkpoint 6 compute heavy Rayleigh integrals. The
% completion time is on the order of several minutes to an hour per
% checkpoint, depending on the specifications of the GPU, possibly even
% longer on a CPU.
%
% All the other checkpoints should complete within a few seconds each.

checkpoint1 = 1;
checkpoint2 = [0 0];
checkpoint3 = [0 0];
checkpoint4 = 0;
checkpoint5 = 0;
checkpoint6 = 0;
checkpoint7 = 0;

% https://www.mathworks.com/help/matlab/matlab_env/
% index-of-code-analyzer-checks.html
%#ok<*UNRCH>

fprintf('==============================================================\n')
fprintf('                           TRANSMIT                           \n')
fprintf('==============================================================\n')

%==========================================================================
% Load hydrohpone data
%==========================================================================
disp('Loading hydrophone data ...')

load(fullfile(PATHS.Data,'ScanData.mat'),'Grid','voltage_mV','medium')

P = voltage_to_pressure(voltage_mV,1/Grid.dt,3,...
    fullfile(PATHS.Data,'sensitivity.mat'));

clear voltage_mV

%==========================================================================
% Compute complex conjugate Fourier
%==========================================================================
disp('Computing conjugate Fourier ...')
P = conj(fft(P,[],3));

%==========================================================================
% Compute angular spectrum
%==========================================================================
disp('Computing angular spectrum ...')

P = angular_spectrum(P,Grid.x,Grid.y);

% Outside the circle k^2 = kx^2 + ky^2 in k-space, all nonzero components
% can be considered to be noise:
f = (0:(Grid.Nt-1))/(Grid.Nt*Grid.dt); % frequency array
k = 2*pi*f/medium.sound_speed; % wavenumber of propagation
P = filter_evanescent_waves(P,Grid,k);
clear f k

%==========================================================================
% Determine scan plane orientation
%==========================================================================
fprintf('\nDetermining scan plane orientation ...\n')

if checkpoint1 == 1
    ScanPlane = determine_scan_plane_orientation(P,Grid,medium);

    save(fullfile(PATHS.Checkpoints,'check1.mat'),'ScanPlane')
elseif checkpoint1 == 2
    load(fullfile(PATHS.Checkpoints,'check1.mat'),'ScanPlane')
else
    fprintf('Pipeline terminated at checkpoint 1.\n')
    return
end

%==========================================================================
% Initial guess scan plane position
%==========================================================================
ScanPlane.x0 = 0;
ScanPlane.y0 = 0;
ScanPlane.z0 = 3e-3; % Aim for a slight underestimation of the value

%==========================================================================
% Convert data back to spatial domain
%==========================================================================

P = angular_spectrum_inverse(P,Grid.x,Grid.y);

for n = 1:N_iterations
    %======================================================================
    % Backward propagation to virtual source plane
    %======================================================================
    fprintf(['\nPropagating backward to virtual source plane ' ...
        '(%.f/%.f) ...\n'],n,N_iterations)

    filename = ['check2-' num2str(n) '.mat'];
    filename = fullfile(PATHS.Checkpoints,filename);

    options.sourceType = 'velocity';
    options.direction  = 'backward';

    if n==1 && showFigures
        options.saveFigure = saveFigures;
    else
        options.saveFigure = false;
    end

    if checkpoint2(n) == 1
        V = rayleigh_propagation(P,Grid,medium,[],ScanPlane,PATHS,options);

        save(filename,'V')
    elseif checkpoint2(n) == 2
        load(filename,'V')
    else
        fprintf(['Pipeline terminated at checkpoint 2 ',...
            '(%.f/%.f).\n'],n,N_iterations)
        return
    end

    % Convert frequency domain data to time domain data
    V = ifft(conj(V),[],3,'symmetric');

    %======================================================================
    % Estimate transducer model parameters
    %======================================================================
    fprintf(['\nEstimating transducer model parameters ',...
        '(%.f/%.f) ...\n'],n,N_iterations)

    filename = ['check3-' num2str(n) '.mat'];
    filename = fullfile(PATHS.Checkpoints,filename);

    if checkpoint3(n) == 1
        [Pmean,ScanPlane,Transducer] = estimate_model_parameters(...
            V,Grid,medium,ScanPlane,PATHS,options);

        save(filename,'ScanPlane','Pmean','Transducer')
    elseif checkpoint3(n) == 2
        load(filename,'ScanPlane','Pmean','Transducer')
    else
        fprintf(['Pipeline terminated at checkpoint 3 ',...
            '(%.f/%.f).\n'],n,N_iterations)
        return
    end

end

% Save data after (re)running the final iteration
if checkpoint2(end) == 1
    save(fullfile(PATHS.Results,'SourceData.mat'),'V')
end
if checkpoint3(end) == 1
    save(fullfile(PATHS.Results,'Transducer.mat'),'Transducer')
    save(fullfile(PATHS.Results,'transformation_scan_plane.mat'),...
        'ScanPlane')
end

%==========================================================================
% Compute transmit impulse response
%==========================================================================
fprintf('\nComputing transmit impulse response ...\n')

options.saveFigure = saveFigures;
options.showFigure = showFigures;

if checkpoint4 == 1
    % Driving pulse
    load(fullfile(PATHS.Data,'TransmitData.mat'),'Transmit');
    y =  Transmit.VoltageSignal;
    Fs = Transmit.SamplingRate;

    [IR,Fs] = compute_transmit_impulse_response(Pmean,1/Grid.dt,y,Fs,...
        medium,ScanPlane,Transducer,PATHS,options);

    save(fullfile(PATHS.Results,'IR_transmit.mat'),'IR','Fs')
elseif checkpoint4 == 0
    fprintf('Pipeline terminated at checkpoint 4.\n')
    return
end

fprintf('\n')
fprintf('==============================================================\n')
fprintf('                           RECEIVE                            \n')
fprintf('==============================================================\n')

%==========================================================================
% Loading receive data from reflection measurement
%==========================================================================
fprintf('\nLoading receive data ...\n')
load(fullfile(PATHS.Data,'ReceiveData.mat'),'Receive')
y  = Receive.RF;
Fs = Receive.SamplingRate;

%==========================================================================
% Determine virtual receiver position and orientation
%==========================================================================
fprintf('\nDetermining virtual receiver position and orientation ...\n')

if checkpoint5 == 1
    VirtualReceiver = determine_virtual_receiver_orientation(y,Fs,...
        medium,Transducer,PATHS,options);

    save(fullfile(PATHS.Results,'transformation_virtual_receiver.mat'),...
        'VirtualReceiver')
elseif checkpoint5 == 2
    load(fullfile(PATHS.Results,'transformation_virtual_receiver.mat'),...
        'VirtualReceiver')
else
    fprintf('Pipeline terminated at checkpoint 5.\n')
    return
end

%==========================================================================
% Forward propagation to virtual sensor plane
%==========================================================================
fprintf('\nPropagating forward to virtual sensor plane ...\n')
options.saveFigure = false;
options.showFigure = showFigures;
options.sourceType = 'pressure';
options.direction  = 'forward';

if checkpoint6 == 1
    [P,rmin] = rayleigh_propagation(...
        P,Grid,medium,ScanPlane,VirtualReceiver,PATHS,options);

    save(fullfile(PATHS.Checkpoints,'check6.mat'),'P','rmin')
elseif checkpoint6 == 2
    load(fullfile(PATHS.Checkpoints,'check6.mat'),'P','rmin')
else
    fprintf('Pipeline terminated at checkpoint 6.\n')
    return
end

% Offset distance used in Rayleigh integral
VirtualReceiver.rmin = rmin;

%==========================================================================
% Revert lens delays and compute average receive pressure
%==========================================================================

% Convert frequency domain data to time domain data
P = ifft(conj(P),[],3,'symmetric');

fprintf(['\nReverting lens delays and '...
    'computing average receive pressure ...\n'])
P = compute_average_receive_pressure(P,Grid,medium,Transducer,options);

%==========================================================================
% Compute transmit impulse response
%==========================================================================

options.saveFigure = saveFigures;

% Make size y Nchannels-by-Nsamples, just like P
y = transpose(y);

fprintf('\nComputing receive impulse response ...\n')
if checkpoint7 == 1
    [IR,Fs] = compute_receive_impulse_response(P,1/Grid.dt,y,Fs,...
        medium,VirtualReceiver,Transducer,PATHS,options);

    save(fullfile(PATHS.Results,'IR_receive.mat'),'IR','Fs')
elseif checkpoint7 == 0
    fprintf('Pipeline terminated at checkpoint 6.\n')
    return
end

fprintf('\n')
fprintf('==============================================================\n')
fprintf('                           COMPLETE                           \n')
fprintf('==============================================================\n')
