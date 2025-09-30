function [IR,Fs] = compute_transmit_impulse_response(p,Fsp,y,Fsy,...
    medium,ScanPlane,Transducer,PATHS,options)
% This function estimates computes the transmit impulse response.
%
% Creates the base plot for Fig. 6d of the arXiv preprint.
%
% input
% - p:   average transmit pressure data (Nsamples-by-1)
% - Fsp: sampling rate of p
% - y:   transmit RF data (Nsamples-by-1)
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

saveFigure = options.saveFigure;
showFigure = options.showFigure;
z0 = ScanPlane.z0;

%==========================================================================
% COMPUTE THE TRANSMIT IMPULSE RESPONSE
%==========================================================================

% Integrated value of the driving pulse [V*s]:
A = sum(y)/Fsy;

% Impulse response in Pa/V/s^2
IR = p/A;

% Estimate the fraction of the transducer surface that is occupied by
% acitve elements
N = Transducer.NumberOfElements;
p = Transducer.Pitch;
w = Transducer.ElementWidth;
fraction = N*w/((N-1)*p + w);

IR = IR/fraction;

% The flat trace before the start of the pulse has been transported to the
% end of the signal (periodic Fourier transform). Remove this part.
M = round(z0/medium.sound_speed*Fsp);
N = length(IR);
IR(N-M+1:N) = [];
N = length(IR);
t1 = (0:(N-1))/Fsp;

% Upsample the signal
Fs = 250e6; % Sampling rate for custom data in PROTEUS
t = t1(1):1/Fs:t1(end);
IR = sinc_interpolation(t1,IR,t);

%==========================================================================
% Visualization
%==========================================================================

if showFigure
    fig = figure;
    plot(t*1e6,IR*1e-9)

    title('Transmit impulse response')
    xlabel('Time (microseconds)')
    ylabel('Pressure (kPa/V/us)')
else
    saveFigure = false;
end

drawnow
if not(saveFigure)
    return
end

%==========================================================================
% Save figure
%==========================================================================

saveFolder = string(fullfile(PATHS.Figures,'fig'));

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

savefig(fig, saveFolder + filesep + "Fig6d.fig")

end
