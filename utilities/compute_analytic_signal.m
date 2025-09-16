function ya = compute_analytic_signal(y)
%COMPUTE_ANALYTIC_SIGNAL - Compute analytic representation
% ya = COMPUTE_ANALTYIC_SIGNAL(y) computes the analytic representation of a
% real-valued signal y.
%
% Output is identical to ya = hilbert(y) from the Signal Processing Toolbox
%
% See https://en.wikipedia.org/wiki/Analytic_signal
%
% See also: hilbert
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

N = length(y);

% Folding frequency at N/2 + 1
M1 = ceil(N/2);      % Frequency index just below folding frequency
M2 = floor(N/2) + 2; % Frequency index just above folding frequency

% Compute analytic signal
Y = fft(y);
Y(M2:N) = 0;
Y(2:M1) = Y(2:M1)*2;
ya = ifft(Y);

end
