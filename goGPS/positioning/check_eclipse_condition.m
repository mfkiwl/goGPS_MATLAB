function [eclipsed] = check_eclipse_condition(time, XS, SP3, sat)

% SYNTAX:
%   [eclipsed] = check_eclipse_condition(time, XS, SP3, sat);
%
% INPUT:
%   time = GPS time
%   XS   = satellite position (X,Y,Z)
%   SP3  = structure containing precise ephemeris data
%   sat  = satellite PRN
%
% OUTPUT:
%   eclipsed = boolean value to define satellite eclipse condition (0: OK, 1: eclipsed)
%
% DESCRIPTION:
%   Check if the input satellite is under eclipse condition.

%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

eclipsed = 0;

t_sun = SP3.t_sun;
X_sun = SP3.X_sun;

[~, q] = min(abs(t_sun - time));
X_sun = X_sun(:,q);

%satellite geocentric position
XS_n = norm(XS);
XS_u = XS / XS_n;

%sun geocentric position
X_sun_n = norm(X_sun);
X_sun_u = X_sun / X_sun_n;

%satellite-sun angle
cosPhi = dot(XS_u, X_sun_u);

%threshold to detect noon/midnight maneuvers
if (~isempty(strfind(SP3.satType{sat},'BLOCK IIA')))
    t = 4.9*pi/180; % maximum yaw rate of 0.098 deg/sec (Kouba, 2009)
elseif (~isempty(strfind(SP3.satType{sat},'BLOCK IIR')))
    t = 2.6*pi/180; % maximum yaw rate of 0.2 deg/sec (Kouba, 2009)
elseif (~isempty(strfind(SP3.satType{sat},'BLOCK IIF')))
    t = 4.35*pi/180; % maximum yaw rate of 0.11 deg/sec (Dilssner, 2010)
end

%shadow crossing affects only BLOCK IIA satellites
shadowCrossing = cosPhi < 0 && XS_n*sqrt(1 - cosPhi^2) < goGNSS.ELL_A_GPS;
if (shadowCrossing && ~isempty(strfind(SP3.satType{sat},'BLOCK IIA')))
    eclipsed = 1;
end

%noon/midnight maneuvers affect all satellites
noonMidnightTurn = acos(abs(cosPhi)) < t;
if (noonMidnightTurn)
    eclipsed = 3;
end
