function P = voltage_to_pressure(V,Fs,DIM,FILENAME)
%VOLTAGE_TO_PRESSURE converts hydrophone voltage data into pressure data.
%
% P = VOLTAGE_TO_PRESSURE(V,Fs,DIM,FILENAME) converts the voltage data V
% (mV) into pressure data P (Pa) by applying the hydrophone sensitivity
% data in FILENAME. V is a multi-dimensional array with dimension DIM
% representing the time axis with sampling rate Fs (Hz).
%
% FILENAME should be a MAT file that contains the variables Frequency (MHz)
% and Sensitivity (mV/MPa), both 1-by-N arrays, where N is the number of
% calibrated hydrophone frequencies.
%
% For Precision Acoustics fibre-optic hydrophones, the file
% read_hydrophone_sensitivity.m can be used to create the variables
% Sensitivity and Frequency.
%
% See also read_hydrophone_sensitivity
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

load(FILENAME,'Frequency','Sensitivity')
fmax = max(Frequency);

% Get the inverse sensitivity:
Sensitivity = 1./Sensitivity;  % [MPa/mV]
Sensitivity = Sensitivity*1e3; % [kPa/mV]

% The polarity of the voltage signals of the hydrophone system are
% inverted. Invert the polarity again:
Sensitivity = -Sensitivity;

% Get the interpolated values of the sensitivity curve
N = size(V,DIM);
f = (0:(N-1))/N*Fs; % Frequency vector [Hz]
fMHz = f/1e6;       % Frequency vector [MHz]
S = interp1(Frequency,Sensitivity,fMHz,'pchip');
S(fMHz>fmax) = 0;

% Reshape the sensitivity array to act on dimension DIM
SHAPE = ones(1,DIM);
SHAPE(DIM) = length(S);
S = reshape(S,SHAPE);

% Apply the sensitivity curve in the frequency domain:
P = fft(V,[],DIM).*S;
P = ifft(P,[],DIM,'symmetric');

% Convert kPa to Pa
P = P*1e3;

end