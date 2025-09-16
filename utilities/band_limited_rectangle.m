function B = band_limited_rectangle(x,xc,dx,W,k)
%BAND_LIMITED_RECTANGLE - Band-limited rectangle function
% Create a band-limited rectangle function with width W, centred at x==xc,
% and band-limited by spatial frequency k.
%
% See equation 6.47 of my thesis
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

N = length(x);

% Fourier transform of a shifted rectangle
kx = 2*pi*(0:(N-1))/(N*dx);
[W,kx] = ndgrid(W,kx);
K = kx.*W/2;
B = W.*sin(K)./K.*exp(-1i*kx.*(xc-x(1)));
B(K==0) = W(K==0);

% Band-limiting
B(kx>k) = 0;

% Transform to spatial domain
B = ifft(B,[],2,'symmetric')/dx;
end