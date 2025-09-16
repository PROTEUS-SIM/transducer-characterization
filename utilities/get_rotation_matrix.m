function R = get_rotation_matrix(theta,axis)
%GET_ROTATION_MATRIX Compute 3D rotation matrix.
% R = GET_ROTATION_MATRIX(theta,axis) computes the 3D rotation matrix of
% theta radians about axis (1, 2, or 3).
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken

if ~isnumeric(theta) || ~isreal(theta) || ~isscalar(theta)
    error('Error: Input argument theta must be a single, real number.')
end

if ~any(axis == [1 2 3])
    error('Error: Input argument axis must have value 1, 2, or 3.')
end

R = eye(3);
R(1:2,1:2) = [cos(theta) -sin(theta); sin(theta) cos(theta)];
permutation = circshift([1 2 3],-axis);
R(permutation,permutation) = R;

end