function phi = get_phase_response_hydrophone(f,method,sensfile)
% Compute the phase response of the hydrophone assuming the hydrophone is a
% minimum-phase system.
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

load(sensfile,'Frequency','Sensitivity')

S1 = griddedInterpolant(Frequency,Sensitivity,method);

fMHz = f/1e6; % Frequency vector [MHz]
phi = compute_minimum_phase(fMHz,S1);
phi = transpose(phi);

end