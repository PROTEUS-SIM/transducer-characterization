function G = compute_spectral_propagator(Nx,dkx,Ny,dky,k,d,...
    averaging,greenType)
%COMPUTE_SPECTRAL_PROPAGATOR computes the spectral propagator G for the
%angular spectrum method.
%
% Nx:  number of grid points in first dimension
% dkx: k-space grid spacing in first dimension
% Ny:  number of grid points in second dimension
% dky: k-space grid spacing in second dimension
% k:   wave number (2*pi/wavelength) [1/m]
% d:   projection distances [m]
% averaging: return averaged spectral propagator (boolean)
% greenType: Green's function type:
% - 'pp': pressure to pressure
% - 'vp': normal velocity to pressure
% - 'pv': pressure to normal velocity
%
% Averaging of the spectral propagator based on:
% Earl G. Williams and J. D. Maynar, Numerical evaluation of the Rayleigh
% integral for planar radiators using the FFT, Journal of the Acoustical
% Society of America, 72(6), 1982, section I.B.3. The case velocity to
% pressure for zero distance is explicitly worked out in the article, the
% other cases are a straightforward extension.
%
% See https://doi.org/10.1121/1.388633
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

% Third dimension corresponds to the propagation axis:
d = reshape(d,1,1,length(d));

% Fourth dimension corresponds to the wave numbers:
k = reshape(k,1,1,1,length(k));

% k-space coordinates:
kx = (ceil(-Nx/2):ceil(Nx/2-1))*dkx;
ky = (ceil(-Ny/2):ceil(Ny/2-1))*dky;
[kx,ky] = ndgrid(kx,ky);

%==========================================================================
% GREEN'S FUNCTIONS EVALUATED AT THE K-SPACE GRID POINTS
%==========================================================================
if ~averaging
    
    % Propagation wavenumber:
    kz = sqrt(k.^2 - kx.^2 - ky.^2);

    switch greenType
        case 'vp' % Normal velocity to pressure
            G = exp(1i*kz.*d).*k./kz;
        case 'pp' % Pressure to pressure
            G = exp(1i*kz.*d);
        case 'pv' % Pressure to normal velocity
            G = kz./k.*exp(1i*kz.*d);
    end
    
    return
end


%==========================================================================
% AVERAGING OF THE GREEN'S FUNCTIONS
%==========================================================================
% Following the integration approach in radial coordinates outlined in
% Williams and Maynar.

% Radial k-space coordinates:
kr = sqrt(kx.^2 + ky.^2);

% Effective radial step size for rectangular grid. Extension of the
% integration approach outlined in section I.B.3 of of Williams and Maynar
% to non-square grids. Approximate the thickness of the integration ring as
% the thickness of a ring defined by the ellipses (kx/dkx)^2 + (ky/dky)^2 =
% K and (kx/dkx + 1)^2 + (ky/dky + 1)^2 = K, where K is a constant.
% (simplifies to dkx and dky along the cartesian axes):
dkr = kr.*((kx/dkx).^2 + (ky/dky).^2).^(-1/2);

% Integration limits for each k-space grid point:
k1 = kr - dkr/2;
k2 = kr + dkr/2;
kz1 = sqrt(k.^2 - k1.^2);
kz2 = sqrt(k.^2 - k2.^2);

switch greenType
    %======================================================================
    % Normal velocity to pressure
    %======================================================================
	case 'vp'
        G = 2./(k2.^2 - k1.^2).*(...
            exp(1i*kz1.*d) - ...
            exp(1i*kz2.*d))./(1i*d).*k;
        
        % Case d == 0
        G(:,:,d==0,:) = 2./(k2.^2 - k1.^2).*(kz1 - kz2).*k...
            .*ones(size(G(:,:,d==0,:)));
        
        % Division by zero at the k-space centre. Replace by the limit kr
        % approaches zero (which is the original function value). In case k
        % == 0, G is zero everywhere, except at the centre. Integration
        % around a small region around the centre yields zero again.
        [i,j] = ind2sub([Nx,Ny],find(kr==0));
        G(i,j,:,:) = exp(1i*abs(k).*d).*sign(k)...
            .*ones(size(G(i,j,:,:)));
                
    case 'pp'
	%======================================================================
    % Pressure to pressure
    %======================================================================
        G = 2./(k2.^2 - k1.^2).*(...
            (1 - 1i*d.*kz1).*exp(1i*kz1.*d) - ...
            (1 - 1i*d.*kz2).*exp(1i*kz2.*d))./d.^2;
        
        % Case d == 0
        G(:,:,d==0,:) = ones(size(G(:,:,d==0,:)));
        
        % Limit kr approaches zero:
        [i,j] = ind2sub([Nx,Ny],find(kr==0));
        G(i,j,:,:) = exp(1i*abs(k).*d)...
            .*ones(size(G(i,j,:,:)));
               
    case 'pv'
	%======================================================================
    % Pressure to normal velocity
    %======================================================================
        G = 2./(k2.^2 - k1.^2).*(...
            (1i*kz2.^2.*d.^2 - 2*kz2.*d - 2i).*exp(1i*kz2.*d) - ...
            (1i*kz1.^2.*d.^2 - 2*kz1.*d - 2i).*exp(1i*kz1.*d))./(k.*d.^3);
        
        % Case d == 0
        G(:,:,d==0,:) = 2./(k2.^2 - k1.^2).*(kz1.^3 - kz2.^3)./(3*k)...
            .*ones(size(G(:,:,d==0,:)));
        
        % Limit kr approaches zero:
        [i,j] = ind2sub([Nx,Ny],find(kr==0));
        G(i,j,:,:) = exp(1i*abs(k).*d)./sign(k)...
            .*ones(size(G(i,j,:,:)));
end

end