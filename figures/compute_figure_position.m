function position = compute_figure_position(figureWidth,figureHeight)
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Get screen properties:
h = get(0);

% Width and height of the figure in pixels:
H = figureHeight*h.ScreenPixelsPerInch;
W = figureWidth *h.ScreenPixelsPerInch;

position = [h.ScreenSize(3)/2-W/2 h.ScreenSize(4)/2-H/2 W H];

end