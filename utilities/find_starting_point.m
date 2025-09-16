function t_start = find_starting_point(...
    y,Fs,upsamplefactor,y_threshold,options)
%FIND_STARTING_POINT - Find the onset time of an ultrasound pulse
%
% T_START = FIND_STARTING_POINT(y,FS,UPSAMPLEFACTOR,Y_THRESHOLD) returns
% the onset time T_START of an ultrasound pulse. T_START is defined as the
% last zero-crossing before the signal exceeds the threshold value
% Y_THRESHOLD.
%
% FS is the sampling rate of the signal and UPSAMPLEFACTOR is an integer to
% upsample the signal for higher precision.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

showFigure = options.showFigure;
saveFigure = options.saveFigure;
saveFolder = options.saveFolder;

% Upsample the signal
Fs2 = Fs*upsamplefactor;
N1 = length(y);
N2 = N1*upsamplefactor;
t1 = (0:(N1-1))/Fs;
t2 = (0:(N2-1))/Fs2;
y = sinc_interpolation_periodic(t1,y,t2);

% Find the first point that exceeds three times the noise level
I = find(y>y_threshold,1);

% Find the last zero-crossing before that point
I = find(y(1:I)<=0,1,'last');

% If no zero-crossing found, make use of assumed periodicity of signal
if isempty(I)
    I = find(y<=0,1,'last');
    Y_start = y(I);
    I = I - N2;
else
    Y_start = y(I);
end

t_start = (I-1)/Fs2;

%--------------------------------------------------------------------------
% Visualization
%--------------------------------------------------------------------------

if showFigure
    displayUnit1 = 1e-6; % (mm)
    displayUnit2 = 1e3;  % (Pa)

    fig = figure;
    plot(t2/displayUnit1,y/displayUnit2)
    hold on
    plot(t_start/displayUnit1,Y_start/displayUnit2,'o')

    title('Pressure at transducer surface')
    xlabel('Time (microseconds)')
    ylabel('Pressure (kPa)')
else
    saveFigure = false;
end

drawnow
if not(saveFigure)
    return
end

%--------------------------------------------------------------------------
% Save figure
%--------------------------------------------------------------------------

if not(isfolder(saveFolder))
    mkdir(saveFolder)
end

savefig(fig, saveFolder + filesep + "Fig6c.fig")

end
