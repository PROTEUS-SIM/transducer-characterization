function [source,sensor] = ...
    rayleigh_integral_frequency(source,sensor,medium)
%RAYLEIGH_INTEGRAL_FREQUENCY computes the frequency domain Rayleigh
%integral for a list of observation points over a set of integration points
%on the transducer surface.
%
% source, a struct with fields:
% - points:      source points [m],                      Nsource-by-Ndim
% - type:        'velocity' or 'pressure'
% - velocities:  normal velocities [m/s], complex,       Nsource-by-Nfreq
% - pressures:   pressure values [Pa], complex,          Nsource-by-Nfreq
% - normal:      normal vector at each source point,     Nsource-by-Ndim
% - weights:     integration element surface area [m^2], Nsource-by-1
% - frequencies: frequency vector [Hz],                  1-by-Nfreq
% - baffle:      'rigid', 'soft', or 'free' (optional)
%
% sensor, a struct with fields:
% - points:      sensor points [m],                      Nsensor-by-Ndim
% - pressures:   pressure values [Pa], complex,          Nsensor-by-Nfreq
% - normal:      normal vector at each sensor point,     Nsensor-by-Ndim
% - weights:     integration element surface area [m^2], Nsensor-by-1
% - offset:      offset distance [m],                    1-by-1
%
% Here, Ndim = 3 is the number of spatial dimensions
%
% medium, a struct with fields:
% - density [kg/m^3]
% - sound_speed [m/s]
% - direction: 'forward' or 'backward'
% - useGPU: set to true for matrix multiplication on the GPU
% 
% RETURNS source, sensor with updated fields:
% sensor.pressures for forward propagation
% source.pressures or source.velocities for backward propagation
%
%
% REQUIRED FIELDS (indicated with X)
%
% propagation: | forward             | backward            |
% ----------------------------------------------------------
% source type: | velocity | pressure | pressure | velocity |
% ==========================================================
% source field |          |          |          |          |
% ----------------------------------------------------------
% points       |    X     |    X     |    X     |    X     |
% type         |    X     |    X     |    X     |    X     | 
% velocities   |    X     |          |          |          |
% pressures    |          |    X     |          |          |
% normal       |          |    X     |          |    X     | 
% weights      |    X     |    X     |          |          |
% frequencies  |    X     |    X     |    X     |    X     |
% baffle       |          |          |          |          | 
% ----------------------------------------------------------             |
% sensor field |          |          |          |          |
% ----------------------------------------------------------
% points       |    X     |    X     |    X     |    X     |
% pressures    |          |          |    X     |    X     |
% normal       |          |          |    X     |    X     |
% weights      |          |          |    X     |    X     |
% offset       |          |          |          |          |
%
% The optional offset distance r0 = sensor.offset defines a translation of
% the origin of the time axis to r0/c0, where c0 is the speed of sound.
% This property can help reduce the required length of the source signals
% in the time domain, thereby reducing computation time and memory.
%
% REFERENCES
%
% Rayleigh integral from Journal of the Acoustical Society of America,
% 72(6), 1982, Eq. 1. Also in: Sapozhnikov et al., J. Acoust. Soc. Am. 138
% (3), 2015. See https://doi.org/10.1121/1.4928396.
% Also in Earl G. Williams and J. D. Maynard, Numerical evaluation of the
% Rayleigh integral for planar radiators using the FFT. See
% https://doi.org/10.1121/1.388633.
%
% See also get_transducer_integration_points, define_grid
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

Nsource = size(source.points,1);      % Number of source points
Nsensor = size(sensor.points,1);      % Number of sensor points
Ndim    = size(source.points,2);      % Number of spatial dimensions (3)
Nfreq   = size(source.frequencies,2); % Number of frequencies

% Matrix dimension corresponding to the spatial coordinates:
DIM = 4;

% Function to cast variables to GPU type (if requested) and reshape:
if medium.useGPU
    reshape_cast = @(A,sz) reshape(gpuArray(A),sz);
else
    reshape_cast = @(A,sz) reshape(A,sz);
end

% Medium properties:
rho0 = medium.density;
c0   = medium.sound_speed;

% Wave number:
k = 2*pi*transpose(source.frequencies)/c0;
k = reshape_cast(k,[1 Nfreq]);

if isfield(sensor,'offset')
    r0 = sensor.offset;
else
    r0 = 0;
end

%==========================================================================
% PROPAGATION KERNEL TYPES
%==========================================================================

% Select integration kernel type:
switch medium.direction
    case 'forward'
        DIM_INT = 3; % Dimension to sum (integrate) over
        switch source.type
            case 'velocity'
                kernel_type = 'fwd_vp';
            case 'pressure'
                kernel_type = 'fwd_pp';
        end
        
    case 'backward'
        DIM_INT = 1; % Dimension to sum (integrate) over
        switch source.type
            case 'velocity'
                kernel_type = 'bwd_pv';
            case 'pressure'
                kernel_type = 'bwd_pp';
        end
end

%==========================================================================
% VECTORS
%==========================================================================

% Compute pairwise distance matrix (Nsensor-by-1-by-Nsource):
r1 = reshape_cast(source.points,[1 1 Nsource Ndim]);
r2 = reshape_cast(sensor.points,[Nsensor 1 1 Ndim]);
r = vecnorm(r1-r2,2,DIM);

% Normal vectors:
if strcmp(kernel_type,'fwd_pp') || strcmp(kernel_type,'bwd_pv')
    
    % Unit vectors pointing from source points to sensor points:
    m12 = (r2 - r1)./r;
    
    % Normal vectors to source:
    n1 = reshape_cast(source.normal,[1 1 Nsource Ndim]);
    n1 = repmat(n1,[Nsensor 1 1 1]);
    
end

if strcmp(kernel_type,'bwd_pp') || strcmp(kernel_type,'bwd_pv')
    
    % Unit vectors pointing from sensor points to source points:
    m21 = (r1 - r2)./r;
    
    % Normal vectors to sensor:
    n2 = reshape(sensor.normal,[Nsensor 1 1 Ndim]);
    n2 = repmat(n2,[1 1 Nsource 1]);
    
end

%==========================================================================
% OBLIQUITY FACTOR
%==========================================================================

% Check baffle type before computing obliquity factor:
if isfield(source,'baffle') && ~strcmp(kernel_type,'fwd_vp')
    warning(['Baffle only supported for forward propagation velocity '...
        'to pressure.'])
end

% Obliquity factor:
if strcmp(kernel_type,'fwd_vp') && ...
        (~isfield(source,'baffle') || strcmp(source.baffle,'rigid'))
    beta = 1;
elseif strcmp(kernel_type,'fwd_vp')
    
    % Cast variables to GPU variables:
    source_points = reshape_cast(source.points, size(source.points));
    sensor_points = reshape_cast(sensor.points, size(sensor.points));
    
    beta = compute_obliquity_factor(source_points, sensor_points, ...
        source.normal, source.baffle);
    beta = reshape_cast(beta,[Nsensor 1 Nsource]);
end

%==========================================================================
% INTEGRATION KERNELS, INTEGRATION VALUES, INTEGRATION WEIGHTS
%==========================================================================

% Forward propagation convolution kernels (Eq. A7 - A10 from Sapozhnikov et
% al.):
switch kernel_type
    case 'fwd_vp'
        S1 = reshape_cast(transpose(source.velocities), [1 Nfreq Nsource]);
        dA = reshape_cast(          source.weights,     [1 1     Nsource]);
        
        K = (-1i*k*rho0*c0)./(2*pi*r).*beta.*exp(1i*k.*(r-r0));
        
    case 'fwd_pp'
        S1 = reshape_cast(transpose(source.pressures),  [1 Nfreq Nsource]);
        dA = reshape_cast(          source.weights,     [1 1     Nsource]);
        
        K = 1/(2*pi)*dot(m12,n1,DIM).*(-1i*k./r + 1./r.^2)...
            .*exp(1i*k.*(r-r0));
        
    case 'bwd_pp'
        S1 = reshape_cast(sensor.pressures,  [Nsensor Nfreq 1]);
        dA = reshape_cast(sensor.weights,    [Nsensor 1     1]);
        
        K = 1/(2*pi)*dot(m21,n2,DIM).*(1i*k./r + 1./r.^2)...
            .*exp(-1i*k.*(r-r0));
        
    case 'bwd_pv'
        S1 = reshape_cast(sensor.pressures,  [Nsensor Nfreq 1]);
        dA = reshape_cast(sensor.weights,    [Nsensor 1     1]);
        
        K = 1./(2i*pi*k*rho0*c0).*(dot(n1,n2,DIM).*(1i*k./r + 1./r.^2) ...
            + dot(m12,n1,DIM).*dot(m21,n2,DIM)...
            .*(3i*k./r + 3./r.^2 - k.^2)).*exp(-1i*k.*(r-r0))./r;
        
end

%==========================================================================
% RAYLEIGH INTEGRAL
%==========================================================================

% Rayleigh integral:
S2 = sum(S1.*K.*dA,DIM_INT);
S2 = gather(S2);

% Assign Rayleigh integral result to the corresponding struct field:
switch kernel_type
    case 'fwd_vp'
        sensor.pressures  = S2; % Nsensor-by-Nfreq
    case 'fwd_pp'
        sensor.pressures  = S2; % Nsensor-by-Nfreq
    case 'bwd_pp'
        source.pressures  = permute(S2,[3 2 1]); % Nsource-by-Nfreq
    case 'bwd_pv'
        source.velocities = permute(S2,[3 2 1]); % Nsource-by-Nfreq
end

end