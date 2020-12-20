function rcm = getMeridianRadiusCurvature(lat_rad)
% get Earth meridian radius of curvature at geodetic latitude
% lat (radians), using GPS ellipsoid (WGS84)
% SYNTAX: 
%    rcm = getMeridianRadiusCurvature(lat_rad)


%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b8
%
%--------------------------------------------------------------------------
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
% Adapted from Octave
%  Written by: Giulio Tagliaferro
%  Contributors:     
%  A list of all the historical goGPS contributors is in CREDITS.nfo
%--------------------------------------------------------------------------
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

a = GPS_SS.ELL_A;
e2 = GPS_SS.ELL_E2;
rcm = a * (1 - e2) / (1 - e2*sin(lat_rad)^2)^(3/2);
end