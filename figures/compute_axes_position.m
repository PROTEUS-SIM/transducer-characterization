function position = compute_axes_position(fig,figureWidth,figureHeight)
% COMPUTE_AXES_POSITION returns the axes position in units of pixels.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Get screen properties:
h = get(0);

% Width and height of the figure in pixels:
H = figureHeight*h.ScreenPixelsPerInch;
W = figureWidth *h.ScreenPixelsPerInch;

position = [fig.Position(3)/2-W/2 fig.Position(4)/2-H/2 W H];

end