function pgrid = define_grid_demos(SimulationParameters,Geometry,varargin)
%DEFINE_GRID_DEMOS - Two-dimensional PROTEUS grid for demos
% This function returns a two-dimensional PROTEUS grid tailored to the
% demos in this repository.
%
% pgrid = DEFINE_GRID_DEMOS(SimulationParameters,Geometry) creates a
% two-dimensional PROTEUS grid by setting the x-dimension of the domain to
% zero.
%
% If Geometry.Axis is specified, the grid is such that the grid lies in the
% plane perpendicular to Geometry.Axis.
%
% pgrid = DEFINE_GRID_DEMOS(SimulationParameters,Geometry,
% expansion_factor) also expands the dimensions of the grid by
% expansion_factor.
%
% pgrid = DEFINE_GRID_DEMOS(SimulationParameters,Geometry,
% expansion_factor, ppwl) also updates the points per wavelength in the
% grid.
%
% See also: define_grid
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

%--------------------------------------------------------------------------
% Input handling
%--------------------------------------------------------------------------

if nargin < 2
    error('Not enough input arguments.')
end

if nargin > 2
    expansion_factor = varargin{1};
else
    expansion_factor = 1;
end

if nargin > 3
    ppwl2 = varargin{2};
else
    ppwl2 = SimulationParameters.PointsPerWavelength;
end

if nargin > 4
    error('Too many input arguments.')
end

%--------------------------------------------------------------------------
% Modify the simulation domain
%--------------------------------------------------------------------------

Geometry.Domain.Xmin = 0;
Geometry.Domain.Xmax = 0;
Geometry.Domain.Ymax = Geometry.Domain.Ymax*expansion_factor;
Geometry.Domain.Ymin = Geometry.Domain.Ymin*expansion_factor;
Geometry.Domain.Zmax = Geometry.Domain.Zmax*expansion_factor;
Geometry.Domain.Zmin = Geometry.Domain.Zmin*expansion_factor;

% Rotate the domain if the requested axis is not the x axis:
if     isfield(Geometry,'Axis') && Geometry.Axis == 2
    Geometry.Domain = permute_domain(Geometry.Domain,[3 1 2]);
elseif isfield(Geometry,'Axis') && Geometry.Axis == 3
    Geometry.Domain = permute_domain(Geometry.Domain,[2 3 1]);
end

%--------------------------------------------------------------------------
% Modify the number of points per wavelength
%--------------------------------------------------------------------------
% Old number of points per wavelength:
ppwl1 = SimulationParameters.PointsPerWavelength;

SimulationParameters.CFL = SimulationParameters.CFL*ppwl1/ppwl2;
SimulationParameters.GridSize = SimulationParameters.GridSize*ppwl1/ppwl2;
SimulationParameters.PointsPerWavelength = ppwl2;

% Except for the demos using k-Wave, no PML or grid optimization is needed
SimulationParameters.MaxPrime = Inf;
SimulationParameters.PMLmin = 0;

%--------------------------------------------------------------------------
% Create the grid
%--------------------------------------------------------------------------

pgrid = define_grid(SimulationParameters,Geometry);

end
