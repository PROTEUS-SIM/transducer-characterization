% This script investigates the factors that may influence the measured
% phase response of the transducer, including the phase response of the
% hydrophone and the phase response of the analog filter. The results are
% described in Section 6.4.2 of my thesis.
%
% Creates the base plots for Fig. 6.10 of my thesis.
%
% Dependencies:
% results/IR_transmit.mat
% results/IR_receive.mat
% data/sensitivity.mat
% data/BLP-15+.txt
%
% Output:
% results/IR_transmit_processed.mat
% results/IR_receive_processed.mat
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

clear
clc
close all

saveFigure = false;

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
addpath(rootDirectory); clear currentDirectory rootDirectory
PATHS = path_setup_calibration();

% Hydrophone sensitivity data file
sensfile = fullfile(PATHS.Data,'sensitivity.mat');

% Analog filter data file
filterfile = fullfile(PATHS.Data,'BLP-15+.txt');

%==========================================================================
% Load data
%==========================================================================

load(fullfile(PATHS.Results,'IR_transmit.mat'),'IR','Fs')
ht = IR;

load(fullfile(PATHS.Results,'IR_receive.mat'),'IR')
hr = IR;
clear IR

Tmax = 2e-6;
N = round(Fs*Tmax);
ht = ht(1:N);
hr = hr(1:N);

ht = ht - mean(ht);
hr = hr - mean(hr);

%==========================================================================
% Compute transfer functions
%==========================================================================

M = 1e4;
Ht = 20*log10(abs(fft(ht,M))/max(abs(fft(ht,M)))); % Transmit
Hr = 20*log10(abs(fft(hr,M))/max(abs(fft(hr,M)))); % Receive

%==========================================================================
% Compensate for the phase of the measurement system
%==========================================================================

f = (0:(N-1))/N*Fs; % Frequency array
interpolationMethod = 'linear';

phi1 = get_phase_response_filter(f,filterfile);
phi2 = get_phase_response_hydrophone(f,interpolationMethod,sensfile);

ht1 = ifft(fft(ht).*exp(-1i*(phi1+phi2)),'symmetric');

phi1 = get_phase_response_filter(f,filterfile);
phi2 = get_phase_response_hydrophone(f,interpolationMethod,sensfile);
 
hr1 = ifft(fft(hr).*exp(1i*(phi1+phi2)),'symmetric');

%==========================================================================
% Show results and save data
%==========================================================================

figure
f = (0:(M-1))/M*Fs; % Frequency array
displayUnit = 1e6; % (Hz)
plot(f/displayUnit,Hr)
hold on
plot(f/displayUnit,Ht)
xlim([0 10])
ylim([-40 5])
xlabel('Frequency (MHz)')
ylabel('Response (dB)')
title('Transfer functions')
legend('Receive','Transmit')
grid on

figure
t = (0:(N-1))/Fs;   % Time array
displayUnit = 1e-6; % (s)
plot(t/displayUnit,hr/max(abs(hilbert(hr))))
hold on
plot(t/displayUnit,ht/max(abs(hilbert(ht))))
ylim([-1 1])
xlabel('Time (microseconds)')
title('Normalized impulse response')
legend('Receive','Transmit')
grid on

figure
plot(t/displayUnit,hr1/max(abs(hilbert(hr1))))
hold on
plot(t/displayUnit,ht1/max(abs(hilbert(ht1))))
ylim([-1 1])
xlabel('Time (microseconds)')
title('Normalized impulse response (phase compensated)')
legend('Receive','Transmit')
grid on

IR = ht1;
save(fullfile(PATHS.Results,'IR_transmit_processed.mat'),'IR','Fs')

IR = hr1;
save(fullfile(PATHS.Results,'IR_receive_processed.mat'),'IR','Fs')
clear IR

%==========================================================================
% Save figures
%==========================================================================

if not(saveFigure)
    return
end

saveFolder = string(fullfile(PATHS.Figures,'fig'));

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

h = get(0,'Children');
savefig(h(3), saveFolder + filesep + "Fig10a.fig")
savefig(h(2), saveFolder + filesep + "Fig10b.fig")
savefig(h(1), saveFolder + filesep + "Fig10c.fig")