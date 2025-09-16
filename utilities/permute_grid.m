function Grid2 = permute_grid(Grid1,permorder)
%PERMUTE_GRID - Permute a PROTEUS grid
% PERMUTE_GRID permutes the coordinates [x, y, z] in GRID1 with
% permutation order PERMORDER.
%
% GRID1 and GRID2 are PROTEUS grids.
% See: https://github.com/PROTEUS-SIM/PROTEUS/blob/main/acoustic-module/define_grid.m
%
% Allowed permutation orders: [1 2 3], [2 3 1], [3 1 2]
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

Grid2 = Grid1;

msg = "Permutation not supported. Allowed values for PERMORDER: " + ...
    "[1 2 3], [2 3 1], [3 1 2]";

if isequal(permorder,[1 2 3])
    return % Identity permutation
elseif isequal(permorder,[2 3 1])
    Grid2.x  = Grid1.y;
    Grid2.y  = Grid1.z;
    Grid2.z  = Grid1.x;
    Grid2.dx = Grid1.dy;
    Grid2.dy = Grid1.dz;
    Grid2.dz = Grid1.dx;
    Grid2.Nx = Grid1.Ny;
    Grid2.Ny = Grid1.Nz;
    Grid2.Nz = Grid1.Nx;
elseif isequal(permorder,[3 1 2])
    Grid2.x  = Grid1.z;
    Grid2.y  = Grid1.x;
    Grid2.z  = Grid1.y;
    Grid2.dx = Grid1.dz;
    Grid2.dy = Grid1.dx;
    Grid2.dz = Grid1.dy;
    Grid2.Nx = Grid1.Nz;
    Grid2.Ny = Grid1.Nx;
    Grid2.Nz = Grid1.Ny;
else
    error(msg)
end

Grid2.full_size = Grid1.full_size(permorder);

if isfield(Grid1,'PML') && isequal(permorder,[2 3 1])
    Grid2.PML.X_SIZE = Grid1.PML.Y_SIZE;
    Grid2.PML.Y_SIZE = Grid1.PML.Z_SIZE;
    Grid2.PML.Z_SIZE = Grid1.PML.X_SIZE;
elseif isfield(Grid1,'PML') && isequal(permorder,[3 1 2])
    Grid2.PML.X_SIZE = Grid1.PML.Z_SIZE;
    Grid2.PML.Y_SIZE = Grid1.PML.X_SIZE;
    Grid2.PML.Z_SIZE = Grid1.PML.Y_SIZE;
end

end