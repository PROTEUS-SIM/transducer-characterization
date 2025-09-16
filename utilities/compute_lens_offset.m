function d_lens = compute_lens_offset(d,Transducer)
%COMPUTE_LENS_OFFSET - compute lens delay in two-way receive data
% d_lens = COMPUTE_LENS_OFFSET(d,Transducer) computes the delay (expressed
% in length units) of an echo pulse onset time caused by the acoustic lens.
%
% d is the round-trip distance from the transducer to the reflector back to
% the transducer.
%
% Transducer, struct with fields:
% - ElevationFocus
% - ElementHeight
%
% See Eq. 6.49 of my thesis.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

F = Transducer.ElevationFocus;
H = Transducer.ElementHeight;

if F<=0
    error('Only positive elevation focus supported')
end

a = H/2;

if d>=2*F
    d_lens = sqrt(F^2+a^2) + sqrt((1-d/F)^2*F^2+a^2) - d;
elseif d>=F
    d_lens = 2*sqrt(F^2+a^2) - 2*F;
elseif isfinite(F)
    d_lens = (1+d/F)*sqrt(F^2+a^2) - sqrt(F^2+(1-d/F)^2*a^2) - d;
else
    d_lens = 0;
end

end
