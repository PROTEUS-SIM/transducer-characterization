function [M,I] = find_reflection_peaks(y,fraction,I0,Npeaks)
%FIND_REFLECTION_PEAKS - Find maxima in signal with reflections
% [M,I] = FIND_REFLECTION_PEAKS(Y,FRACTION,I0,NPEAKS) finds the values M
% and the indices I of the maxima in a signal Y containing multiple
% internal reflections from a reflector: i.e., a decaying pulse train.
%
% First, the envelope of the signal is computed. Next, the first maximum
% from index I0 is searched. Next, the first point after which the signal
% decays to FRACTION is searched. Next, the following maximum is searched.
% This is repeated NPEAKS times.
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Compute signal envelope
A = abs(compute_analytic_signal(y));
N = length(A);

M = zeros(1,Npeaks);
I = zeros(1,Npeaks);

for k = 1:Npeaks

    % Find first maximum from index I0:
    [M1,I1] = max(A(I0:end)); I1 = I0+I1-1;
    
    % Find the first point where the envelope decays to fraction of the
    % peak value:
    I0 = find(A(I1:N)<A(I1)*fraction,1); I0 = I1+I0-1;

    if not(isempty(I1))
        I(k) = I1;
        M(k) = M1;
    end

end

end
