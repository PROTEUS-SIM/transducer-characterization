function P = angular_spectrum_inverse(A,x,y)
%ANGULAR_SPECTRUM_INVERSE Compute the inverse two-dimensional Fourier
%transform along the first two dimensions.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

xs = x([find(x>=0) find(x<0)]);
ys = y([find(y>=0) find(y<0)]);

A = ifftshift(A,1);
A = ifftshift(A,2);
P = ifft2(A);
P = P([find(xs<0) find(xs>=0)],[find(ys<0) find(ys>=0)],:,:);

end