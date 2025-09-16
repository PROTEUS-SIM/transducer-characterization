function VirtualReceiver = determine_virtual_receiver_orientation(...
    y,Fs,medium,Transducer,PATHS,options)
% This function computes the orientation and position of the virtual
% receiver. This corresponds to the following step in the blue part of Fig.
% 6.4 of my thesis:
% - Estimate virtual receiver position and orientation
%
% Creates the base plots for panels b and c of Fig. 6.7 of my thesis.
%
% input
% - y:  Receive RF data (Nt-by-Nchannels)
% - Fs: Receive RF sampling rate
%
% output
% - VirtualReceiver, struct with fields thetaX, thetaY, thetaZ, x0, y0, z0
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%==========================================================================
% Parameters
%==========================================================================

% Lower bounds for the start of the reflected signal
t1 = 40e-6; % Lower bound 1 (start of flat noisy trace)
t2 = 60e-6; % Lower bound 2 (just before signal onset)

fc = 0.1e6; % Cut-off frequency for high-pass filter

% As the reflector and the transducer are not perfectly aligned, the
% incident wave is not plane over the full range of transducer elements.
% Between the following elements the received wave can be considered plane
% (see plot of receive amplitude vs element number):
element_range = 10:70;

%==========================================================================
% Show data
%==========================================================================

saveFigure = options.saveFigure;
showFigure = options.showFigure;

% Time vector:
N = size(y,1);
t = (0:(N-1))/Fs;

if showFigure
    figure
    plot(max(y(t>t2,:)))
    title('Receive amplitude')
    xlabel('Element number')
    ylabel('Maximum amplitude')

    fig2 = figure;
    imagesc(t(t>t2)*1e6,1:Transducer.NumberOfElements,y(t>t2,:)')
    title('Receive data')
    xlabel('Receive time (microseconds)')
    ylabel('Element number')

    hold on
end

% Find the angle between transducer and incoming plane wave
[thetaY, delays] = find_angle_receive(t(t>t2),y(t>t2,:),...
    medium.sound_speed,Transducer,element_range,options);

%==========================================================================
% Compensate for the angle in the received plane wave
%==========================================================================

% Shift channel signals (by modulation in the frequency domain) such that
% all channel signals line up:
f = transpose((0:(N-1))/N*Fs);
y = fft(y);
y = y.*exp(2i*pi*f.*delays);
y = ifft(y,'symmetric');

% Shifting the channels is circular and moves data from the beginning of
% the signal to the end. Compute the latest time that is free from this
% data:
t3 = t(end)-max(delays);

% Average signals to obtain a better signal-to-noise ratio:
y = mean(y(:,element_range),2);

% The MHz-range ultrasound signal is superposed on a kHz-range electronic
% signal. Even though this signal has an amplitude two to three orders of
% magnitude lower than the ultrasound signal, it makes it difficult to
% determine the starting point of the ultrasound signal. Filter out the
% lower-frequency signal with a Butterworth highpass filter.
[b,a] = butter(3,fc/(Fs/2),'high'); % Filter coefficients:
y = filter(b,a,y); % Zero-phase filter

if showFigure
    fig3 = figure;
    plot(t(t>t1&t<=t3)*1e6,y(t>t1&t<=t3))
    title('Average receive signal (angle compensated)')
    xlabel('Receive time (microseconds)')
    ylabel('Amplitude')
    hold on
end

%==========================================================================
% Determine the starting point of the first reflection
%==========================================================================

% Get an onset threshold based on the noise level
th = std(y(t>t1&t<t2))*5;

% Find the last zero crossing before crossing the threshold value
I = find(transpose(y)>th&t>t2,1);
I = find(y(1:I)<=0,1,'last');
t0 = t(I);

if showFigure
    plot(t0*1e6,y(I),'o')
else
    saveFigure = false;
end

% Compute the delay of the wave front caused by the lens
D_lens = compute_lens_offset(t0*medium.sound_speed,Transducer);

% Centre of the virtual sensor
z0 = t0*medium.sound_speed - D_lens;
x0 = z0*sind(thetaY)/(1+cosd(thetaY));

VirtualReceiver.thetaX = 0;
VirtualReceiver.thetaY = thetaY;
VirtualReceiver.thetaZ = 0;

VirtualReceiver.x0 = x0;
VirtualReceiver.y0 = 0;
VirtualReceiver.z0 = z0;
VirtualReceiver.t0 = t0;

drawnow
if not(saveFigure)
    return
end

%==========================================================================
% Save figures
%==========================================================================

saveFolder = string(fullfile(PATHS.Figures,'fig'));

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

savefig(fig2, saveFolder + filesep + "Fig7b.fig")
savefig(fig3, saveFolder + filesep + "Fig7c.fig")

end