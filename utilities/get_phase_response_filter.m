function phi = get_phase_response_filter(f,filterfile)
%GET_FILTER_PHASE_RESPONSE returns the phase response of an analog low pass
%filter.
% PHI = GET_FILTER_PHASE_RESPONSE(F) returns the phase PHI at the query
% frequencies F.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

fc = 15e6; % Cutoff frequency of the filter (Hz)

[~,~,~,~, FrequencyMHz2, GroupDelaynsec] = read_filter_data(filterfile);

% Find the phase through integration of the group delay
fcMHz = fc/1e6;    % Cutoff frequency in MHz
dfMHz = fcMHz/1e3; % Step size for integration

fMHz = (0:dfMHz:FrequencyMHz2(end));
tau  = interp1(FrequencyMHz2,GroupDelaynsec/1e3,fMHz,'pchip');
phi  = -cumtrapz(tau)*dfMHz*2*pi;

% Get the phase at the query frequencies
phi = griddedInterpolant(fMHz,phi);
phi = phi(f/1e6);

end
