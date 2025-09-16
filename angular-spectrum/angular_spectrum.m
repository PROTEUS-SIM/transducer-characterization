function A = angular_spectrum(P,x,y)
%ANGULAR_SPECTRUM Compute the two-dimensional Fourier transform along the
%first two dimensions.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

P = P([find(x>=0) find(x<0)],[find(y>=0) find(y<0)],:,:);
A = fft2(P);
A = fftshift(A,1);
A = fftshift(A,2);

end