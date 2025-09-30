function G = filter_evanescent_waves(G,Grid,varargin)
%FILTER_EVANESCENT_WAVES - Filter evanescent wave components
% This function filters evanescent wave components from angular spectrum
% data.
%
% If the angular spectrum data G has size (NX,NY,NK), use
% G = FILTER_EVANESCENT_WAVES(G,GRID,K), where K are the wave vectors for
% the filter (row array).
%
% If the angular spectrum data G has size (NX,NY,ND,NK), use
% G = FILTER_EVANESCENT_WAVES(G,GRID,K,D), where K are the wave vectors for
% the filter (row array) and D is are the distances of each plane (row
% array).
%
% This file is part of the transducer-characterization project, licensed
% under the GNU Lesser General Public License v3.0 (LGPL-3.0).
% See the LICENSE file for further details.
% Copyright (C) 2025 Nathan Blanken


Grid = get_kspace(Grid);

if nargin==2
    error('Not enough input arguments.');
elseif nargin==3
    k = varargin{1};
    [kx,ky,k] = ndgrid(Grid.kx,Grid.ky,k);
elseif nargin==4
    k = varargin{1};
    d = varargin{2};
    [kx,ky,~,k] = ndgrid(Grid.kx,Grid.ky,d,k);
else
    error('Too many input arguments.')
end

dk = max(Grid.dkx,Grid.dky); % Margin for filter
G(kx.^2 + ky.^2 >= (k-dk).^2) = 0;
G(k==0) = 0; % k = 0 component

end
