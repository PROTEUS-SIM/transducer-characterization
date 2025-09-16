function Rp = compute_reflection_coefficient(p,medium_attenuating)
%COMPUTE_REFLECTION_COEFFICIENT - compute the pressure reflection
%coefficient
% Rp = COMPUTE_REFLECTION_COEFFICIENT(P,medium_attenuating) computes the
% pressure reflection coefficient of a reflecting plate from an array of
% pulse amplitudes P based on a 1D model of pulse decay.
%
% The pulse amplitudes P correspond to the multiple reflections from the
% reflecting plate.
%
% Setting medium_attenuating = true will return a more general solution,
% but requires three pulse amplitudes (1-by-3 array).
%
% Setting medium_attenuating = false will return a less general solution,
% but only requires two pulse amplitudes (1-by-2 array) and is therefore
% less sensitive to experimental error.
%
% 1D model for decay of reflected pulses. p0 is the incident pulse. p1, p2,
% p3, etc. are pulses reflected back to the transducer.
%
% p1 = Rp*p0
% p2 = -Rp^1*TI*exp(-2*a*D)*p0
% p3 = -Rp^3*TI*exp(-4*a*D)*p0
% pn = -Rp^(2*n-3)*TI*(exp-2*(n-1)*a*D)*p0 for n > 1
% 
% Rp: pressure reflection coefficient
% TI: intensity transmission coefficient
% D:  plate thickness (m)
% a:  medium attenuation (Np/m)
%
% Rp = (x-1)/(x+1)
% TI = 4*x/(x+1)^2
% with x = Z1/Z2 the ratio between first and second acoustic impedance
%
% Ratios between pulse amplitudes:
% abs(p3/p2) = Rp^2*exp(-2*a*D)
% abs(p2/p1) = TI*exp(-2*a*D)
% abs(p3/p2)/abs(p2/p1) = Rp^2/TI = (x-1)^2/(4*x)
% With A = abs(p3/p2)/abs(p2/p1)
% x = 2*A + 1 +/- 2*sqrt(A^2 + A)
%
% If medium not attenuation (a = 0):
% abs(p2/p1) = TI = 4*x/(x+1)^2
% With A = abs(p1/p2)
% x = 2*A - 1 +/- 2*sqrt(A^2 - A)
% 
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Compute ratio between acoustic impedances x = Z2/Z1
if medium_attenuating
    A = abs(p(3)/p(2))/abs(p(2)/p(1));
    x = 1 + 2*A + 2*sqrt(A^2 + A);
else
    A = abs(p(1)/p(2));
    x = 2*A - 1 + 2*sqrt(A^2 - A);
end

% Compute reflection coefficient:
Rp = (x-1)/(x+1);

end
