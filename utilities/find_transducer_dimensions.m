function [W,H,xc,yc] = find_transducer_dimensions(P,Parameters)
%FIND_TRANSDUCER_DIMENSIONS - Extract the width, height and centre of the
%transducer surface.
%
% [W,H,xc,yc] = FIND_TRANSDUCER_DIMENSIONS extracts the width W, height H
% and centre (xc,yc) of the transducer surface based on a single-frequency
% pressure/velocity plane P (Nx-by-Ny).
%
% Parameters should contain the following fields:
% - x:   Lateral coordinates
% - y:   Elevation coordinates
% - dx:  Lateral grid spacing
% - dy:  Elevation grid spacing
% - W:   Initial guess for width
% - H:   Initial guess for height
% - xc:  Initial guess for centre
% - yc:  Initial guess for centre
% - sound_speed: (m/s)
% - frequency:   (Hz)
% - showFit: boolean
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

disp('Fitting rectangular aperture in lateral direction ...')
Y = mean(abs(P),2);
Y = transpose(Y);
Y = Y/max(Y);

Parameters.Axis = 1;
[a,W,xc,E] = find_fit(Y,Parameters);

fprintf('RMS error:         %5.2f %%\n',  E/a*100)
fprintf('Transducer width:  %5.2f mm\n',  W*1e3)
fprintf('Transducer centre: %5.2f mm\n\n',xc*1e3)

disp('Fitting rectangular aperture in elevation direction ...')
Y = mean(abs(P),1);
Y = Y/max(Y);

Parameters.Axis = 2;
[a,H,yc,E] = find_fit(Y,Parameters);

fprintf('RMS error:         %5.2f %%\n',  E/a*100)
fprintf('Transducer width:  %5.2f mm\n',  H*1e3)
fprintf('Transducer centre: %5.2f mm\n\n',yc*1e3)

end

function [a,W,xc,e] = find_fit(Y,Parameters)

% Toggle axis
if Parameters.Axis == 1
    x   = Parameters.x;
    dx  = Parameters.dx;
    W0  = Parameters.W;  % Initial guess for width
    xc0 = Parameters.xc; % Initial guess for centre
    saveName = "Fig6a1.fig";
elseif Parameters.Axis == 2
    x   = Parameters.y;
    dx  = Parameters.dy;
    W0  = Parameters.H;  % Initial guess for width
    xc0 = Parameters.yc; % Initial guess for centre
    saveName = "Fig6a2.fig";
else
    error('Invalid axis')
end

% Currently the speed of sound in water is used to compute the maximum
% spatial frequency. However, using a higher speed of sound results in even
% better fits (speed of sound in PZT perhaps?)
factor = 1;
c0 = Parameters.sound_speed*factor;
f0 = Parameters.frequency;

% Maximum spatial frequency
kmax = 2*pi*f0/c0;

% Function handle for the band-limited rectangle functions with the
% rectangle centre xc and width W as fitting parameters:
B = @(xc,W) band_limited_rectangle(x,xc,dx,W,kmax);

% The optimal amplitude A for a given xc, and W can be expressed
% analytically: A = sum(Y*B)/sum(B*B) (L2 linear regression):
A = @(xc,W) sum(Y.*B(xc,W),2)./sum(B(xc,W).*B(xc,W),2);

% Fitting function (amplitude multiplied by shape):
F = @(xc,W) A(xc,W).*B(xc,W);

% We have taken the absolute of A. Therefore it may seem fair to take the
% absolute of the fitting curve B as well. However, this will make the
% optimization problem nonconvex.

% Error (L2 objective function):
E = @(p) sum((F(p(1),p(2)) - Y).^2);

% Optimization:
[p,e] = fminsearch(E, [xc0 W0]);
xc = p(1);              % Rectangle centre
W  = p(2);              % Rectangle width
a  = A(p(1),p(2));      % Amplitude of the best fit
e  = sqrt(e/length(Y)); % RMS error

%--------------------------------------------------------------------------
% Visualization
%--------------------------------------------------------------------------

if Parameters.showFit
    displayUnit = 1e-3; % (m)
    fig = figure;
    plot(x/displayUnit,Y/a)
    hold on
    plot(x/displayUnit,F(xc,W)/a)

    if Parameters.Axis == 1
        title('Best fit to aperture (Lateral)')
        xlabel('x (mm)')
    else
        title('Best fit to aperture (Elevation)')
        xlabel('y (mm)')
    end
else
    Parameters.saveFit = false;
end

if not(Parameters.saveFit)
    return
end

%--------------------------------------------------------------------------
% Save figure
%--------------------------------------------------------------------------

saveFolder = Parameters.saveFolder;

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

savefig(fig, saveFolder + filesep + saveName)

end