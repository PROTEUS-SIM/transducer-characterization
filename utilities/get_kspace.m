function Grid = get_kspace(Grid)
%GET_KSPACE - Add k-space to Grid struct
% GRID = GET_KSPACE(GRID) computes the kspace vectors of the struct GRID
% and adds these vectors to the struct.
%
% The k-space vectors are defined such that the zero component corresponds
% to element ceil((N+1)/2), in line with fftshift.
%
% See also: fftshift
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

Grid.kx  = 2*pi*(ceil(-Grid.Nx/2):ceil(Grid.Nx/2-1))/(Grid.Nx*Grid.dx);
Grid.ky  = 2*pi*(ceil(-Grid.Ny/2):ceil(Grid.Ny/2-1))/(Grid.Ny*Grid.dy);
Grid.kz  = 2*pi*(ceil(-Grid.Nz/2):ceil(Grid.Nz/2-1))/(Grid.Nz*Grid.dz);
Grid.dkx = 2*pi/(Grid.Nx*Grid.dx);
Grid.dky = 2*pi/(Grid.Ny*Grid.dy);
Grid.dkz = 2*pi/(Grid.Nz*Grid.dz);

end