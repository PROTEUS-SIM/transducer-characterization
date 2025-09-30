function [F,yc] = find_elevation_focus(P,Parameters)
%FIND_ELEVATION_FOCUS - Find elevation focus based on lens delays.
%
% Returns elevation focus F and lens centre line yc.
%
% Pfft is temporal frequency domain pressure or velocity data
% (Nx-by-Ny-by-Nfrequencies).
%
% See Eq. 48 of the arXiv preprint
%
% Parameters, struct with fields:
% frequencies: full frequency array
% frequency: current frequency component (used for plotting)
% sound_speed: (m/s)
% x:    lateral coordinates
% y:  elevation coordinates
% x1:   left transducer boundary
% x2:  right transducer boundary
% y1: bottom transducer boundary
% y2:    top transducer boundary
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

f = Parameters.frequencies;

% Frequency component used for plotting:
[~,I0] = min(abs(f-Parameters.frequency));

% Fit a parabolic curve for each frequency component
M  = length(f);
p1 = zeros(1,M); % First polynomial coefficients
p2 = zeros(1,M); % Second polynomial coefficients
A  = zeros(1,M); % Weights based on amplitude

for n = 1:length(f)

    Parameters.frequency = f(n);

    if n == I0 && Parameters.showFit
        Parameters.plotting = true;
    else
        Parameters.plotting  = false;
    end

    [p1(n), p2(n), A(n)] = fit_polynomial_delays(P(:,:,n),Parameters);

end

% Take a weighted average of the polynomial fit
A = A/sum(A);
p1 = sum(A.*p1);
p2 = sum(A.*p2);

c0 = Parameters.sound_speed;
F  = -1/(2*p1*c0); % Elevation focus (m)
yc = -p2/(2*p1);  % Lens centre line (m)

end

function [p1, p2, A] = fit_polynomial_delays(P,Parameters)

f0 = Parameters.frequency;
c0 = Parameters.sound_speed;

% Due to the band-limited nature of the data, the delays near the edges of
% the transducer are not a reliable representation of the lens delays.
% Select a subdomain sufficiently far from the edges.
margin = c0/f0;
x1 = Parameters.x1 + margin; [~,ix1] = min(abs(Parameters.x - x1));
x2 = Parameters.x2 - margin; [~,ix2] = min(abs(Parameters.x - x2));
y1 = Parameters.y1 + margin; [~,iy1] = min(abs(Parameters.y - y1));
y2 = Parameters.y2 - margin; [~,iy2] = min(abs(Parameters.y - y2));

x = Parameters.x(ix1:ix2);
y = Parameters.y(iy1:iy2);

% Subdomain y-coordinates and signal values
P = P(ix1:ix2,iy1:iy2);
[~,Y] = ndgrid(x,y);

% Compute lens delays for each point
phi = angle(P);
phi = unwrap(phi,[],2);
delay = -phi/(2*pi*f0);

% Weigth for currenty frequency component
A = sum(abs(P(:)));

% Make sure there are at least three distinct y-values to ensure a properly
% conditioned fit.
if (iy2-iy1+1)>=3 && (ix2-ix1+1)>=1
    % Parabolic fit:
    pfit = polyfit(Y(:),delay(:),2);
    p1 = pfit(1);
    p2 = pfit(2);

else
    % If fit is not properly conditioned:
    p1 = 0;
    p2 = 0;
end

%--------------------------------------------------------------------------
% Visualization
%--------------------------------------------------------------------------

if Parameters.plotting
    fig1 = figure;
    imagesc(x*1e3,y*1e3,delay')
    title('Delays (cropped to aperture)')
    xlabel('Lateral (mm)')
    ylabel('Elevation (mm)')

    fig2 = figure;
    plot(Y(:)*1e3,delay(:)*1e6,'.')
    hold on
    plot(y*1e3,polyval(pfit,y)*1e6)
    hold off
    title('Best fit to delays')
    xlabel('Elevation (mm)')
    ylabel('Delay (microseconds)')
else
    Parameters.saveFit = false;
end

if not(Parameters.saveFit)
    return
end

%--------------------------------------------------------------------------
% Save figures
%--------------------------------------------------------------------------

saveFolder = Parameters.saveFolder;

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

savefig(fig1, saveFolder + filesep + "Fig6b.fig")
savefig(fig2, saveFolder + filesep + "Fig6b1.fig")

end
