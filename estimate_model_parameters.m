function [P,ScanPlane,Transducer] = estimate_model_parameters(...
    Y,Grid,medium,ScanPlane,PATHS,options)
% This function estimates transducer model parameters. This function
% corresponds to the following parts of Fig. 6.4 of my thesis:
% - Estimate (x0, y0); Estimate transducer dimensions W and H
% - Estimate lens focus F
% - Revert lens delays
% - Compute average transmit velocity
% - Updated x0, y0, z0
% - Determine z0
% - Dimensions W, H; Focus F
%
% Creates the base plots for Fig. 6.6a-c of my thesis.
%
% Input:
% - Y: frequency-domain pressure or velocity data (depending on
%      options.sourceType) (Nx-by-Ny-by-Nsamples)
% - ScanPlane, struct with fields thetaX, thetaY, thetaZ, x0, y0, z0
%
% Output:
% - P: average pressure over transducer surface (Nsamples-by-1)
% - ScanPlane, struct with fields thetaX, thetaY, thetaZ, x0, y0, z0
% - Transducer, struct with transducer properties, PROTEUS style
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%==========================================================================
% Parameters
%==========================================================================

% Frequency range (Hz) for finding the elevation focus One row per search
% range Scan between 0 and 5 MHz, but do not include the frequencies
% between 3.2 and 3.8 MHz to filter out the spurious mode of oscillation
% (see Section 6.3.4 of my thesis).
freqRange = [0.1 3.2; 3.8 5]*1e6;

displayUnit = 1e-3; % (m)

sourceType = options.sourceType;
saveFigure = options.saveFigure;

saveFolder = string(fullfile(PATHS.Figures,'fig'));
options.saveFolder = saveFolder;

%==========================================================================
% LOAD PRESSURE/VELOCITY DATA
%==========================================================================

if strcmp(sourceType,'pressure')
    P = Y;
elseif strcmp(sourceType,'velocity')
    % Convert velocity data to plane-wave equivalent pressure
    P = Y*medium.density*medium.sound_speed;
end

Fs = 1/Grid.dt;
f  = (0:(Grid.Nt-1))/Grid.Nt*Fs;

Pfft = fft(P,[],3)/Fs;

%==========================================================================
% FIND TRANSDUCER DIMENSIONS AND CENTRE
%==========================================================================

% Load manufacturer details for initial estimates:
load(fullfile(PATHS.Data,'TransmitData.mat'),'Trans');
p = Trans.spacingMm;
w = Trans.elementWidth*Trans.spacingMm/Trans.spacing;
N = Trans.numelements;
Parameters.W  = (N*p + w)*1e-3; % Transducer width
Parameters.H  = Trans.elevationApertureMm*1e-3; % Transducer height
Parameters.xc = 0; % Transducer centre
Parameters.yc = 0; % Transducer centre

% Centre frequency and corresponding frequency array index
Parameters.frequency = Trans.frequency*1e6;
[~,I] = min(abs(f-Parameters.frequency));

Parameters.sound_speed = medium.sound_speed;
Parameters.x = Grid.x; Parameters.dx = Grid.dx;
Parameters.y = Grid.y; Parameters.dy = Grid.dy;

Parameters.showFit    = options.showFigure;
Parameters.saveFit    = options.saveFigure;
Parameters.saveFolder = options.saveFolder;
[W,H,xc,yc] = find_transducer_dimensions(Pfft(:,:,I),Parameters);

x1 = xc-W/2; Parameters.x1 = x1;
x2 = xc+W/2; Parameters.x2 = x2;
y1 = yc-H/2; Parameters.y1 = y1;
y2 = yc+H/2; Parameters.y2 = y2;

if options.showFigure
    fig = figure;
    imagesc(Grid.x/displayUnit,Grid.y/displayUnit,abs(Pfft(:,:,I))')
    hold on
    plot([x1 x1 x2 x2 x1]/displayUnit,[y1 y2 y2 y1 y1]/displayUnit,'r')
    title('Best fit to aperture (boundary)')
    xlabel('Lateral coordinate (mm)')
    ylabel('Elevation coordinate (mm)')
else
    saveFigure = false;
end

%==========================================================================
% ESTIMATE GEOMETRICAL ELEVATION FOCAL DISTANCE
%==========================================================================

% Frequency components to be included in the scan:
I = [];
for n = 1:size(freqRange,1)
    I = [I find(f>freqRange(n,1)&f<freqRange(n,2))]; %#ok
end

Parameters.frequencies = f(I);

F = find_elevation_focus(Pfft(:,:,I),Parameters);

disp(['The geometrical elevation focus is at ' num2str(F) ' m.' newline])

%==========================================================================
% REVERSE LENS EFFECT
%==========================================================================

% Undo lens delays
c0    = medium.sound_speed;
[~,y] = ndgrid(Grid.x,Grid.y);
y     = y - yc;
tau   = -y.^2/(2*F*c0) + (H/2)^2/(2*F*c0);

f    = reshape(f,1,1,length(f));
Pfft = Pfft.*exp(2i*pi*f.*tau);
P    = ifft(Pfft,[],3,'symmetric')*Fs;

%==========================================================================
% AVERAGE THE SIGNALS OVER THE SURFACE OF THE TRANSDUCER
%==========================================================================

[~,ix1] = min(abs(Grid.x-(xc-W/2)));
[~,ix2] = min(abs(Grid.x-(xc+W/2)));
[~,iy1] = min(abs(Grid.y-(yc-H/2)));
[~,iy2] = min(abs(Grid.y-(yc+H/2)));

P = squeeze(mean(P(ix1:ix2,iy1:iy2,:),[1,2]));

% The flat trace before the start of the pulse has been transported to the
% end of the signal (periodic Fourier transform). Use this part to estimate
% the noise level to set a threshold value.
N = length(P);
M = round(ScanPlane.z0/medium.sound_speed*Fs);
P_threshold = std(P((N-M+1):N))*5;

% Find the onset time of the signal
t0 = find_starting_point(P,Fs,10,P_threshold,options);

% Update scan plane centre coordinates:
ScanPlane.x0 = ScanPlane.x0 - xc;
ScanPlane.y0 = ScanPlane.y0 - yc;
ScanPlane.z0 = ScanPlane.z0 + c0*t0;

disp(['The distance to the transducer surface is ' ...
    num2str(ScanPlane.z0) ' m.'])

%==========================================================================
% Transducer object for PROTEUS
%==========================================================================

% TRANSMIT DATA
% Load transmit data
load(fullfile(PATHS.Data,'TransmitData.mat'),'Trans');

% Get kerf from Verasonics data
p = Trans.spacingMm;
w = Trans.elementWidth*Trans.spacingMm/Trans.spacing;
N = Trans.numelements;
kerf = (p-w)*1e-3; % (m)

% Update element width and pitch based on measured transducer width
w = (W-(N-1)*kerf)/N;
p = (W-w)/(N-1);

% Create Transducer object for PROTEUS
Transducer.NumberOfElements = N;
Transducer.Pitch            = p;
Transducer.ElementWidth     = w;
Transducer.ElementHeight    = H;
Transducer.ElevationFocus   = F;

drawnow
if not(saveFigure)
    return
end

%==========================================================================
% Save figure
%==========================================================================

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

savefig(fig, saveFolder + filesep + "Fig6a.fig")

end
