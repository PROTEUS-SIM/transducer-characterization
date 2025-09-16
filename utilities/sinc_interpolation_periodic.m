function Vq = sinc_interpolation_periodic(x,V,xq)
%SINC_INTERPOLATION - 1-D periodic sinc interpolation
% Performs sinc interpolation to find Vq, the values of the band-limited 
% function V=F(X) at the query points Xq. This function assumes the signal
% V is periodic, i.e.: V(N+1) = V, where N is the length of the signal.
%
% Based on:
% T. Schanze, "Sinc interpolation of discrete periodic signals,"
% IEEE Transactions on Signal Processing, vol. 43, no. 6, pp. 1502-1503,
% June 1995, doi: 10.1109/78.388863.
% https://www.doi.org/10.1109/78.388863
%
% Wise, E. S., Cox, B. T., Jaros, J., & Treeby, B. E. (2019). Representing 
% arbitrary acoustic source and sensor distributions in Fourier collocation
% methods. The Journal of the Acoustical Society of America, 146(1), 
% 278-288.
% https://doi.org/10.1121/1.5116132
%
% See also: sinc_interpolation, evaluate_delta_function (part of the
% PROTEUS repository)
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

dx = mean(diff(x));
N = length(x);

[xq,x] = ndgrid(xq,x);

if mod(N,2)==0
    % Even grid size (Eq. 12):
    b = sin(pi*(x-xq)./dx)./(N.*tan(pi*(x-xq)./(N.*dx)));
else
    % Odd grid size (Eq. 10):
    b  = sin(pi*(x-xq)./dx)./(N.*sin(pi*(x-xq)./(N.*dx)));
end

% Define case where x = xq (0/0):
b(x==xq) = 1;

if isrow(V)
    V = V';
    wasrow = true;
else
    wasrow = false;
end

Vq = b * V;

if wasrow
    Vq = Vq';
end

end
