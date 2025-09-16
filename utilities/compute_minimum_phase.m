function phi = compute_minimum_phase(omega,G)
% Compute the minimum phase of a frequency response G according to John
% Bechhoefer; KramersKronig, Bode, and the meaning of zero. Am. J. Phys. 1
% October 2011; 79 (10): 10531059. https://doi.org/10.1119/1.3614039
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Integration variable
dnu    = 1e-3;
nu_min = -4;
nu_max = 4;
nu     = nu_min:dnu:nu_max;

omega = reshape(omega,length(omega),1);

omega1 = omega*exp(nu);
G = G(omega1);
M = log(abs(G));

dM = zeros(size(M));
dM(:,2:end)   = 1/2*diff(M,[],2);
dM(:,1:end-1) = dM(:,1:end-1) + 1/2*diff(M,[],2);

% Convolution kernel
f  = 2/pi^2*log(coth(abs(nu)/2));
f(f==Inf)=0;

phi = -pi/2*sum(f.*dM,2);

end