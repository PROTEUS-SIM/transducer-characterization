function [IR,Fs] = compute_receive_impulse_response(p,Fsp,y,Fsy,...
    medium,VirtualReceiver,Transducer,PATHS,options)
% This function computes the receive impulse response of the transducer.
% This corresponds to the following step in the blue part of Fig. 6.4 of my
% thesis:
% - Compute receive impulse response
% The results are described in Section 6.3.2 of my thesis.
%
% Creates the base plots for panels d, e, and f of Fig. 6.7 of my thesis.
%
% input
% - p:   average receive pressure data (Nchannels-by-Nsamples)
% - Fsp: sampling rate of p
% - y:   receive RF data (Nchannels-by-Nsamples)
% - Fsy: sampling rate of y
%
% output
% - IR: receive impulse response (1-by-Nsamples)
% - Fs: sampling rate of IR
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%==========================================================================
% Parameters
%==========================================================================

% Between these elements the received wave can be considered plane (see
% also determine_virtual_receiver_orientation.m)
element_range = 10:70;

% Approximate centre frequency of transducer
f0 = 2.5e6;

% Noise floor estimates
noise_level_y = -70; % dB
noise_level_p = -50; % dB

%==========================================================================
% Load and show RF data
%==========================================================================

% Set up time array pressure data
N = size(p,2);
t1 = (0:(N-1))/Fsp + VirtualReceiver.rmin/medium.sound_speed;

% Set up time array receive data
N = size(y,2);
t = (0:(N-1))/Fsy;

nElements = size(y,1);

% Find shared time domain and interpolate pressure signal
T1 = max(min(t),min(t1));
T2 = min(max(t),max(t1));
y = y(:,t>T1&t<T2); t = t(t>T1&t<T2);
p = sinc_interpolation(t1,transpose(p),t);

figure
imagesc(t*1e6,1:nElements,y)
title('Receive data (voltage)')
xlabel('Receive time (microseconds)')
ylabel('Element number')

fig2 = figure;
imagesc(t*1e6,1:nElements,p)
title('Receive data (pressure)')
xlabel('Time (microseconds)')
ylabel('Element number')

%==========================================================================
% Compensate for transmit voltage difference
%==========================================================================
% The hydrophone data and the transducer receive data were not obtained
% with the same transmit voltage. Scale the pressure data to match the
% experimental setting of the receive data.

% Get the voltage used for the pressure measurments
load(fullfile(PATHS.Data,'TransmitData.mat'),'Transmit');
voltage1 = max(Transmit.VoltageSignal);

% Get the voltage used for the reflected signal measurements
load(fullfile(PATHS.Data,'ReceiveData.mat'),'Transmit')
voltage2 = max(Transmit.VoltageSignal);

p = p*voltage2/voltage1;

%==========================================================================
% Compensate for angle
%==========================================================================

thetaY = VirtualReceiver.thetaY;
t0 = VirtualReceiver.t0;

% Transducer element positions
x_el = (1:Transducer.NumberOfElements)*Transducer.Pitch;
x_el = x_el - mean(x_el);

% Set up frequency array
N = length(t);
f = (0:(N-1))/N*Fsy;

% Align all RF element signals and average
delays = -tand(thetaY)*x_el/medium.sound_speed;
delays = delays(:);

y = fft(y,[],2); y = y.*exp(2i*pi*delays.*f);
p = fft(p,[],2); p = p.*exp(2i*pi*delays.*f);

y = ifft(y,[],2,'symmetric');
p = ifft(p,[],2,'symmetric');

y = mean(y(element_range,:));
p = mean(p(element_range,:));

y = y - mean(y);
p = p - mean(p);

%==========================================================================
% Find period between reflections and reflection coefficient
%==========================================================================

% The receive signal contains multiple pulses corresponding to the
% reflection from the reflector's front surface and multiple reflections
% from the reflector's back surface.

I0 = round((t0-t(1))*Fsy) + 1; % Start index of first pulse
A = abs(compute_analytic_signal(y)); % Envelope
[M,I] = find_reflection_peaks(y,0.2,I0,5);
T = mean(diff(I))/Fsy; % Period between reflections
Rp = compute_reflection_coefficient(M,false); % Reflection coefficient

figure
plot(t*1e6,y)
hold on
plot(t*1e6,A)
plot(t(I)*1e6,M,'o')
title('Reflection peaks')
xlabel('Receive time (microseconds)')

% In the simulated receive pressure data, we assumed perfect reflection.
% Compensate for the reflection coefficient:
p = Rp*p;

% Early truncation
% Truncating the signal before deconvolution results in cleaner plots, but
% comes at the expense of a slight loss of information and causality.
if isfield(options,'earlyTruncate') && options.earlyTruncate == true
    % Truncate signal with 3-cycle tapered window
    M = floor(T*Fsy) + I0; % Start of the second pulse
    A = zeros(1,N);
    A(1:M) = tukeywin(M,3/(T*f0));
    y = A.*y;
    plot(t*1e6,y)
end

%==========================================================================
% Deconvolve in the frequency domain
%==========================================================================

Y = fft(y); YdB = 20*log10(abs(Y)./max(abs(Y)));
P = fft(p); PdB = 20*log10(abs(P)./max(abs(P)));

fig4 = figure;
plot(f/1e6,YdB)
hold on
plot(f/1e6,PdB)
title('Deconvolution in frequency domain')
xlabel('Frequency (MHz)')
ylabel('Amplitude spectrum (dB)')
legend('Voltage signal', 'Pressure signal')

% Compute the SNR for both signals
noise_level_y = 10^(noise_level_y/20)*max(abs(Y));
noise_level_p = 10^(noise_level_p/20)*max(abs(P));

SNR_y = (abs(Y)/noise_level_y).^2;
SNR_p = (abs(P)/noise_level_p).^2;

% Modified Wiener deconvolution. See Eq. 6.52 of my thesis.
H = Y./P./(1 + 1./SNR_y + 1./SNR_p);

HdB = 20*log10(abs(H)./max(abs(H)));
plot(f/1e6,HdB)
legend('Voltage signal', 'Pressure signal','Transfer function')

% Convert transfer function to impulse response
IR = ifft(H,'symmetric')*Fsy;

%==========================================================================
% Truncate and upsample impulse response for use in PROTEUS
%==========================================================================

% Upsample the signal
t1 = t-t(1);
Fs = 250e6; % Sampling rate for custom data in PROTEUS
t = 0:1/Fs:t1(end);
IR = sinc_interpolation(t1,IR,t);

% Remove the part of the signal that contains additional reflections
IR1 = IR; IR = IR(t<T);
t1  =  t;  t =  t(t<T);

fig5 = figure;
plot(t1*1e6,IR1*1e-6)
hold on
plot(t*1e6,IR*1e-6)
title('Receive impulse response')
xlabel('Time (microseconds)')
ylabel('Response (ADCL/Pa/us)')

drawnow
if not(options.saveFigure)
    return
end

%==========================================================================
% Save figures
%==========================================================================

saveFolder = string(fullfile(PATHS.Figures,'fig'));

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

savefig(fig2, saveFolder + filesep + "Fig7d.fig")
savefig(fig4, saveFolder + filesep + "Fig7e.fig")
savefig(fig5, saveFolder + filesep + "Fig7f.fig")

end
