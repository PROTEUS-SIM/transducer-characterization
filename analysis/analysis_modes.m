% This script allows to investigate the spurious mode of oscillation. The
% results are described in Section 6.3.4 of my thesis.
%
% Creates the base plots for Fig. 6.9 of my thesis.
%
% Dependencies:
% results/Transducer
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%#ok<*UNRCH>

clear
clc
close all

saveFigure = false;

averagingType = 'none';

pauseTimeSeconds = 0;
showVideo = true;
saveVideo = false;
savename = 'TransducerSurface.avi';
displayUnit = 1e-3; % (m)

if saveFigure == true
    showVideo = false;
    saveVideo = false;
    averagingType = 'lateral';
end

i1 = 117;
i2 = 123;

t1 = [1.2 5 11]*1e-6;
t2 = t1 + 0.15e-6;

% Add root directory to path and run path setup:
currentDirectory = fileparts(mfilename('fullpath'));
rootDirectory = fileparts(currentDirectory);
addpath(rootDirectory); clear currentDirectory rootDirectory
PATHS = path_setup_characterization();

%==========================================================================
% LOAD PRESSURE DATA AND REVERSE LENS DELAYS
%==========================================================================
filename = fullfile(PATHS.Data,'ScanData.mat');
load(filename,'Grid','medium')
load(fullfile(PATHS.Results,'SourceData.mat'),'V')
Fs = 1/Grid.dt;
f  = (0:(Grid.Nt-1))/Grid.Nt*Fs;
Vfft = fft(V,[],3);

% Dimensions and centre of transducer
load(fullfile(PATHS.Results,'Transducer.mat'),'Transducer')

F  = Transducer.ElevationFocus;
H  = Transducer.ElementHeight;
W  = Transducer.Pitch*(Transducer.NumberOfElements-1) ...
    + Transducer.ElementWidth;

% Undo lens delays
c0    = medium.sound_speed;
[~,y] = ndgrid(Grid.x,Grid.y);
tau   = -y.^2/(2*F*c0) + (H/2)^2/(2*F*c0);

f    = reshape(f,1,1,length(f));
Vfft = Vfft.*exp(2i*pi*f.*tau);
V    = ifft(Vfft,[],3,'symmetric');

N = size(V,3);
t = (0:(N-1))/Fs;

figure
plot(t*1e6,squeeze(V(i1,i2,:))*1e3)
xlabel('Time (microseconds)')
ylabel('Velocity (mm/s)')


%==========================================================================
% SHOW A VIDEO OF THE OSCILLATING TRANSDUCER SURFACE
%==========================================================================

Vmax = max(abs(V(:)));

if showVideo
    Nmax = size(V,3);
else
    Nmax = 0;
end

if saveVideo && not(isfolder('videos'))
    mkdir('videos')
end

if saveVideo
    v = VideoWriter(fullfile('videos',savename));
    open(v)
end

h = figure();

% Static plot
if showVideo == false
    
    % Transducer boundaries
    [~,ix1] = min(abs(Grid.x + W/2));
    [~,ix2] = min(abs(Grid.x - W/2));
    
    imagesc(Grid.y*1e3,t*1e6,squeeze(mean(V(ix1:ix2,:,:),1))')
    xlabel('z (mm)')
    ylabel('Time (microseconds)')
    
    figure
    
    % Indices of selected time points
    [~,I1] = min(abs(t-t1'),[],2);
    [~,I2] = min(abs(t-t2'),[],2);
    y = Grid.y;
    for n = 1:length(t1)
        subplot(3,1,n)
        plot(Grid.y*1e3,mean(V(ix1:ix2,:,I1(n)),1)*1e3,'r')
        hold on
        plot(Grid.y*1e3,mean(V(ix1:ix2,:,I2(n)),1)*1e3,'k--')
        title([num2str(t1(n)*1e6) ' microseconds'])
        
        xlabel('y (mm)')
        ylabel('Velocity (mm/s)')
    end
    
end

for k = 1:Nmax
    
    switch averagingType
        case 'none'
            imagesc(Grid.x/displayUnit, Grid.y/displayUnit, ...
                squeeze(V(:,:,k))')
            clim([-Vmax Vmax]*1.5)
            axis equal
            xlabel('Lateral (mm)')
            ylabel('Elevation (mm)')
        case 'lateral'
            plot(Grid.y/displayUnit,mean(V(:,:,k),1))
            ylim([-Vmax Vmax]*1.5)
            xlabel('Elevation (mm)')
            ylabel('Velocity (m/s)')
        case 'elevation'
            plot(Grid.x/displayUnit,mean(V(:,:,k),2))
            ylim([-Vmax Vmax]*1.5)
            xlabel('Lateral (mm)')
            ylabel('Velocity (m/s)')
    end
    
    drawnow
    
    if saveVideo == true
        frame = getframe(h);
        writeVideo(v,frame);
    end
    
    pause(pauseTimeSeconds)
    
end

if saveVideo == true
    close(v)
end

%==========================================================================
% SHOW DISPERSION RELATION
%==========================================================================

% Disregard the section of the section containing thickness thickness
% expander mode.
f0 = 2.5e6; % Centre frequency
Nc = 5;     % Approximate number of cycles of the short pulse
i_start = round(Nc*Fs/f0)+1;
show_dispersion_curve(V(:,:,i_start:end),averagingType,Grid)


if not(saveFigure)
    return
end

saveFolder = string(fullfile(PATHS.Figures,'fig'));

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

h = get(0,'Children');
savefig(h(4), saveFolder + filesep + "Fig9a.fig")
savefig(h(3), saveFolder + filesep + "Fig9b.fig")
savefig(h(2), saveFolder + filesep + "Fig9c.fig")
savefig(h(1), saveFolder + filesep + "Fig9d.fig")


function show_dispersion_curve(V,averagingType,Grid)

switch averagingType
    case 'none'
        % Default to averaging in the lateral direction
        V  = squeeze(mean(V,1));
        dk = Grid.dky;
    case 'lateral'
        V  = squeeze(mean(V,1));
        dk = Grid.dky;
    case 'elevation'
        V  = squeeze(mean(V,2));
        dk = Grid.dkx;
end

Fs = 1/Grid.dt;

% Compute both the spatial and the temporal Fourier transform
Nk = 1e3;
N  = 1e3;
V = fft(V,N,2);
V = fft(V,Nk,1);

f = (0:(N-1))*Fs/N;
w = 2*pi*f;
k = (0:(Nk-1))*dk;

figure
h = axes();
imagesc(k/1e6,squeeze(w)/1e6,abs(V)');

set(h,'YDir','normal')

xlabel('k (\mum^{-1})')
ylabel('\omega_0 (MHz)')
title('Dispersion relation')

xlim([0 N*dk/2]/1e6)
ylim([0 2*pi*Fs/2]/1e6)

end

