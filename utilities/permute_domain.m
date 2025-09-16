function Domain2 = permute_domain(Domain1,permorder)
%PERMUTE_DOMAIN - permute a PROTEUS Geometry domain
%
% DOMAIN1 and DOMAIN2 are PROTEUS domains. A PROTEUS domain is a field of a
% PROTEUS Geometry struct.
% See: https://github.com/PROTEUS-SIM/PROTEUS/blob/main/GUIfunctions/compute_simulation_domain.m
%
% Allowed permutation orders: [1 2 3], [2 3 1], [3 1 2]
%
% This file is part of the transducer-calibration project, licensed under
% the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

Domain2 = Domain1;

msg = "Permutation not supported. Allowed values for PERMORDER: " + ...
    "[1 2 3], [2 3 1], [3 1 2]";

if isequal(permorder,[1 2 3])
    return % Identity permutation
elseif isequal(permorder,[2 3 1])
    Domain2.Xmax  = Domain1.Ymax;
    Domain2.Ymax  = Domain1.Zmax;
    Domain2.Zmax  = Domain1.Xmax;
    Domain2.Xmin  = Domain1.Ymin;
    Domain2.Ymin  = Domain1.Zmin;
    Domain2.Zmin  = Domain1.Xmin;
elseif isequal(permorder,[3 1 2])
    Domain2.Xmax  = Domain1.Zmax;
    Domain2.Ymax  = Domain1.Xmax;
    Domain2.Zmax  = Domain1.Ymax;
    Domain2.Xmin  = Domain1.Zmin;
    Domain2.Ymin  = Domain1.Xmin;
    Domain2.Zmin  = Domain1.Ymin;
else
    error(msg)
end

Domain2.TransducerSurface    = Domain1.TransducerSurface(:,   permorder);
Domain2.TransducerProjection = Domain1.TransducerProjection(:,permorder);
Domain2.Vertices             = Domain1.Vertices(:,            permorder);

end