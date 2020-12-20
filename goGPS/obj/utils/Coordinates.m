%   CLASS Coordinates
% =========================================================================
%
% DESCRIPTION
%   Class to manage coordinates
%
% EXAMPLE
%   pos = Coordinates();
%
% CONSTRUCTOR SYNTAX
%   t = Coordinates.fromXYZ(XYZ);
%
% FOR A LIST OF CONSTANTs and METHODS use doc Coordinates
%
% COMMENTS
% The class stores arrays of time, not just a single element,
% it has been designed this way because MATLAB works best on arrays
%

%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b8
%
%--------------------------------------------------------------------------
%  Copyright (C) 2009-2019 Mirko Reguzzoni, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti, ...
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


classdef Coordinates < Exportable & handle
    
    properties (Constant, GetAccess = public)
        ELL_A = GPS_SS.ELL_A;       % Ellipsoid semi-major axis
        ELL_E = GPS_SS.ELL_E;       % Ellipsoid eccentricity
        E2 = Coordinates.ELL_E^2;   % Square of the ellipsoidal eccentricity
        
        DEG2RAD = pi/180;           % Convert degree to radius
        RAD2DEG = 180/pi;           % Convert radius to degree
    end
    
    properties (SetAccess = public, GetAccess = public) % set permission have been changed from private to public (Giulio)
        name = '';                  % Name of the point (not yet used extensively)
        time = GPS_Time             % Position time
        xyz = []                    % Coordinates are stored in meters in as cartesian XYZ ECEF [m]
        v_xyz = []                  % Coordinates velocities XYZ ECEF  [m / year]
        precision = 0.0001          % 3D limit [m] to check the equivalence among coordinates
        Cxx = [] 
    end
        
    % =========================================================================
    %    
    % =========================================================================
    
    methods
        function this = Coordinates()
            % Constructor
            %
            % SYNTAX
            %   t = Coordinates();            
        end
                        
        function copyFrom(this, pos)
            % Copy from an object of the same type
            %
            % SYNTAX
            %   this.copyFrom(pos)
            this.xyz = pos.xyz;
            this.Cxx = pos.Cxx;
            this.time = pos.time.getCopy;
            try % legacy support
                this.v_xyz = pos.v_xyz;
            catch
                this.v_xyz = [];
            end
        end
        
        function copy = getCopy(this)
            % Get a copy of this
            %
            % SYNTAX
            %   copy = getCopy(this)
            copy = Coordinates();
            copy.copyFrom(this);
        end
        
        function this = append(this, pos)
            % Append a Coordinates object into the this
            %
            % SYNTAX
            %   this = append(this, pos)
            
            this.xyz = [this.xyz; pos.xyz];
            if ~isempty(this.Cxx) &&  ~isempty(pos.Cxx)
                this.Cxx = cat(this.Cxx, pos.Cxx, 3);                
            end
            this.time.append(pos.time);
        end
        
        function rem(this, idx)
            % Remove coordinates into the this
            %
            % SYNTAX
            %   this = rem(this, idx)
            this.xyz(idx,:) = [];
            if ~isempty(this.Cxx) 
                this.Cxx(:,:,idx) = [];
            end
            this.time.remEpoch(idx);
        end
        
        function check(this)
            % Raise an error and empty the coordinates if time and
            % coordinates have different dimensions
            if size(this.xyz) ~= this.time.length
                fprintf('WARNING Coordinates with unconsistent dimension found\n');
                this.time = GPS_Time;
                this.xyz = [];
                this.Cxx = [];
            end
        end
    end
        
    % =========================================================================
    %    GETTERS
    % =========================================================================
    
    methods
        function time = getTime(this)
            % Get the time of the coordinates
            %
            % SYNTAX
            %   time = this.getTime()
            
            time = this.time.getCopy;
        end

        function coo = getMedianPos(sta_list)
            % get the median of the coordinates
            %
            % SYNTAX
            %   coo = getMedianPos(this)
            coo = Coordinates();
            for c = 1 : numel(sta_list)
                this = sta_list(c);
                try
                    if isempty(this.time) || this.time.isEmpty
                        coo(c) = Coordinates.fromXYZ(median(this.xyz, 1, 'omitnan'));
                    else
                        coo(c) = Coordinates.fromXYZ(median(this.xyz, 1, 'omitnan'), this.time.getCentralTime);
                    end
                catch
                    Core.getLogger.addWarning('No data found in coordinate object');
                    coo(c) = Coordinates();
                end
            end
        end
        
        function coo = getElement(this, id_el)
            % get a copy of the coordinates relative to a specific id (or epoch)
            %
            % UNFINISHED FUNCITON: 
            %   covariances and precisions are not extracted
            %
            % SYNTAX
            %   coo = getElement(this, id_el)
            try
                if isempty(this.time) || this.time.isEmpty
                    coo = Coordinates.fromXYZ(this.xyz(id_el, :));
                else
                    coo = Coordinates.fromXYZ(this.xyz(id_el, :), this.time.getEpoch(id_el));
                end
            catch
                Core.getLogger.addWarning('Coordinates get Element out of bound');
                coo = Coordinates();
            end           
        end
                
        function [xyz, y, z, time] = getXYZ(this)
            % Get Coordinates as cartesian Earth Centered Earth Fixed Coordinates
            %
            % OUTPUT
            %   xyz      = coordinates     [m]
            %
            % SYNTAX 
            %   [xyz, time] = this.getXYZ()
            %   [x, y, z, time] = this.getXYZ()
           
            if nargout >= 3
                xyz = this.xyz(:, 1);
                y =   this.xyz(:, 2);
                z =   this.xyz(:, 3);
                if isempty(this.time)
                    this.time = GPS_Time();
                else
                    time = this.time.getCopy;
                end
            else
                xyz = this.xyz;
                if isempty(this.time)
                    this.time = GPS_Time();
                else
                    y = this.time.getCopy;
                end
            end
        end
        
        function [east, north, utm_zone, time] = getUTM(this)
            % Get Coordinates as UTM coordinates
            %
            % OUTPUT
            %   east     = east Coordinates  [m]
            %   north    = north Coordinates [m]
            %   utm_zone = UTM zone          [4char]
            %
            % SYNTAX 
            %   [east, north, utm_zone, time] = this.getUTM();
            
            [lat, lon] = this.getGeodetic();
            [east, north, utm_zone] = this.geod2plan(lat, lon);
            if nargout == 4
                time = this.time.getCopy;
            end
        end
                       
        function [east, north, up, utm_zone, time] = getENU(this, theta)
            % Get Coordinates as UTM ENUs coordinates
            %
            % INPUT
            %   theta   planar rotation angle [degree -360:360]
            %
            % OUTPUT
            %   east     = east Coordinates  [m]
            %   north    = north Coordinates [m]
            %   up       = up                [m]
            %   utm_zone = UTM zone          [4char]
            %
            % SYNTAX 
            %   [east, north, up, utm_zone, time] = this.getENU(<theta>);
            %   [enu, utm_zone, time] = this.getENU(<theta>);
            
            [lat, lon, up] = this.getGeodetic();
            [east, north, utm_zone] = this.geod2plan(lat, lon);
            
            if nargin == 2 && ~isempty(theta)
                tmp = [east(:) north(:)] * [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
                east = tmp(:,1);
                north = tmp(:,2)';
            end            
            if nargout <= 3
                east = [east(:) north(:) up(:)];
                north = utm_zone;
                if nargout == 3
                    up = this.time.getCopy;
                end
            end
            if nargout == 5
                time = this.time.getCopy;
            end
        end
        
        function [lat, lon, h_ellips, h_ortho] = getGeodetic(this)
            % Get Coordinates as Geodetic coordinates
            %
            % OUTPUT
            %   lat      = latitude                      [rad]
            %   lon      = longitude                     [rad]
            %   h_ellips = ellipsoidal height            [m]
            %   h_ortho  = orthometric height            [m]
            %
            % SYNTAX 
            %   [lat, lon, h_ellips, h_ortho] = this.getGeodetic()
           
            if nargout > 2
                [lat, lon, h_ellips] = Coordinates.cart2geod(this.xyz);
                if nargout == 4
                    h_ortho = h_ellips - this.getOrthometricCorrFromLatLon(lat, lon);
                end
            else
                [lat, lon] = Coordinates.cart2geod(this.xyz);
                if nargout <= 1
                    lat = [lat, lon];
                end
            end
        end
        
        function [lat_geoc, lon, h] = getGeocentric(this)
            % Get Coordinates as Geodetic coordinates
            %
            % OUTPUT
            %   east     = east Coordinates [m]
            %   north    = north Coordinates [m]
            %   utm_zone = UTM zone
            %
            % SYNTAX 
            %   [lat, lon, h] = this.getGeocentric();
           
            [lat_geoc, lon, h] = Coordinates.cart2geoc(this.xyz);
        end
        
        function ondu = getOrthometricCorrection(this)
            % Get Orthometric correction from the geoid loaded in Core
            %
            % OUTPUT
            %   ondu     = geoid ondulation [m]
            %
            % SYNTAX 
            %   ondu = getOrthometricCorrection(this)
            
            [lat, lon] = this.getGeodetic();
            ondu = this.getOrthometricCorrFromLatLon(lat, lon);
        end
        
        function [loc_enu, v_enu, id_ok] = getLocal(this, ref_pos)
            % Get Coordinates as Local coordinates with respect to ref_pos
            %
            % OUTPUT
            %   loc     = xyz local coordinates
            %
            % SYNTAX 
            %   loc = this.getLocal(ref_pos)
           
            xyz_ref = ref_pos.getXYZ;
            [xyz_this, time] = this.getXYZ;
            baseline = xyz_this - repmat(xyz_ref, size(xyz_this,1),1);
            [loc_enu, rot_mat] = Coordinates.cart2loca(xyz_ref, baseline);
            if nargout > 1
                if size(xyz_this, 1) > 3
                    if isempty(this.v_xyz)
                        v_enu = [0 0 0];
                        
                        for c = 1:3
                            [~, id_ok(:,c), trend] = strongFilterStaticData(loc_enu(:,c), 0.8, 7);
                            v_enu(c) = (trend(end) - trend(1)) / (time.last.getMatlabTime - time.first.getMatlabTime) * 365; % m / year
                        end
                        this.v_xyz = v_enu * rot_mat;
                    else
                        v_xyz = this.v_xyz; %#ok<PROPLC>
                        v_enu = v_xyz * rot_mat'; %#ok<PROPLC>
                    end
                else
                    id_ok = true(size(xyz_ref,1));
                    v_enu = [nan nan nan];
                end
                
            end
        end
        
        function status = isEmpty(this)
            % Return the status of emptyness of the object
            %
            % SYNTAX
            %   status this.isEmpty();
            status = isempty(this.xyz);
        end
        
        function n_el = length(this)
            % Return the number of coordinates present in the object
            %
            % SYNTAX
            %   n_el = this.length();
            n_el = size(this.xyz, 1);
        end
        
        function dist = ellDistanceTo(this, coo)
            % return the distance on the ellipsoid between the object
            % coordinates and the coordinate coo using vincenty's formula
            % (https://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf)
            %
            % SYNTAX
            %     dist = this.distanceTo(coo)
            
            % reduced latitude
            a = this.ELL_A;
            b = sqrt((1 - this.E2)*a^2) ;
            f = (a - b)/a;
            [lat1, lon1] = this.getGeodetic();
            [lat2, lon2] = coo.getGeodetic();
            U1 = atan((1-f)*tan(lat1));
            U2 = atan((1-f)*tan(lat2));
            L = lon2 - lon1;
            lambda_old = 10000;
            lambda= L;
            
            while abs(lambda - lambda_old) > 1e-12
                sin_sigma = sqrt((cos(U2) * sin(lambda))^2 + (cos(U1)*sin(U2) - sin(U1)*cos(U2)*cos(lambda))^2);
                cos_sigma = sin(U1)*sin(U2) + cos(U1)*cos(U2)*cos(lambda);
                sigma = atan2(sin_sigma,cos_sigma);
                sin_alpha = cos(U1)*cos(U2)*sin(lambda)/sin_sigma;
                cos2_alpha = 1 - sin_alpha^2;
                cos_2sigmam = cos_sigma - 2*sin(U1)*sin(U2)/cos2_alpha;
                C = f/16*cos2_alpha*(4 + f*(4 -3*cos2_alpha));
                lambda_old = lambda;
                lambda = L + (1 - C)*f*sin_alpha*(sigma + C*sin_sigma*(cos_2sigmam + C*cos_sigma*(-1 + 2 * cos_2sigmam^2)));
            end
            u2 = cos2_alpha*(a^2 - b^2)/b^2;
            A = 1 + u2/16384*(4096 + u2*(-768 + u2*(320-175*u2)));
            B = u2/1024*(256 +u2*(-128 +u2*(74-47*u2)));
            Dsigma = B*sin_sigma*(cos_2sigmam + 1/4*B*(cos_sigma*(-1 + 2*cos_2sigmam^2) -B/6*cos_2sigmam*(-3 + 4 * sin_sigma^2)*(- 3 + 4 * cos_2sigmam^2)));
            dist = b*A*(sigma -Dsigma);            
        end
    end
    
    % =========================================================================
    %    SETTERS
    % =========================================================================
    
    methods
        function setTime(this, time)
            % Set the time of the coordinates
            %
            % SYNTAX
            %   this.setTime(time)
            
            if time.length > 3
                this.time = time.getNominalTime(time.getRate/2);
            else
                this.time = time.getCopy;
            end
            if this.time.length > size(this.xyz, 1)
                this.time.getEpoch(1 : size(this.xyz, 1));
                fprintf('The set coordinates time is larger than the number of positions stored\nCutting time\nDebug from Coordinates.setTime()');
            elseif this.time.length > size(this.xyz, 1)
                this.xyz = this.xyz(1 : this.time.length, :);
                this.time.getEpoch(1 : size(this.xyz, 1));
                fprintf('The set coordinates time is smaller than the number of positions stored\Cutting positions\nDebug from Coordinates.setTime()');                
            end

        end
        
        function setPosXYZ(this, xyz, y, z)
            % Set the Coordinates
            %
            % SYNTAX
            %   this.setPosXYZ(xyz)
            %   this.setPosXYZ(x, y, z)
            
            if nargin == 4
                xyz = [xyz(:) y(:) z(:)];
            end
            this.xyz = xyz;
        end
        
        function empty(this)
            % Empty the object
            %
            % SYNTAX
            %   this.empty();
            this.xyz = [];
            this.time = GPS_Time();
        end
        
        function addOffset(this, enu_offset, time_start, time_stop)
            % Add the offset (of the anntenna) to a set of coordinates
            %
            % SYNTAX
            %   this.addOffset(enu_offset, <time_start>, <time_stop>)
            if nargin == 4
                id_ok = (this.time >= time_start & this.time >= time_stop);
            else
                id_ok = 1 : size(this.xyz,1);
            end
            
            xyz_offset = Coordinates.local2cart(this.getElement(id_ok).getMedianPos.xyz, enu_offset);
            this.xyz = this.xyz + repmat(xyz_offset, size(this.xyz,1), 1);
        end
    end
    
    % =========================================================================
    %    STATIC CONSTRUCTOR
    % =========================================================================
    methods (Access = 'public', Static)
        function this = fromXYZ(xyz, y, z, time)
            % Set the Coordinates from XYZ coordinates
            %
            % SYNTAX
            %   this = Coordinates.fromXYZ(xyz)
            %   this = Coordinates.fromXYZ(x, y, z)
            
            this = Coordinates;
            if nargin > 2
                xyz = [xyz(:) y(:) z(:)];
                if nargin == 3
                    time = GPS_Time();
                end
            else
                if nargin == 2
                    time = y;
                else
                    time = GPS_Time();
                end
            end
            
            this.setPosXYZ(xyz);
            this.setTime(time);
        end
        
        function this = fromStringXYZ(xyz_string, time)
            % Set the Coordinates from XYZ coordinates (String)
            %
            % SYNTAX
            %   this = Coordinates.fromStringXYZ(xyz)
            
            this = Coordinates;
            xyz = sscanf(xyz_string, '%f%f%f')';
            if numel(xyz) ~= 3
                xyz = [0 0 0];
            end
            if nargin == 1
                time = GPS_Time();
            end
            this.setPosXYZ(xyz);
            this.setTime(time);
        end
        
        function this = fromGeodetic(lat, lon, h_ellips, h_ortho, time)
            % Set the Coordinates from Geodetic coordinates
            %
            % INPUT
            %   lat, lon [rad]
            %
            % SYNTAX
            %   this = Coordinates.fromGeodetic(phi, lam, h_ellips);
            %   this = Coordinates.fromGeodetic(phi, lam, [], h_ortho);
            
            this = Coordinates;
            if nargin >= 4
                h_ellips = h_ortho + this.getOrthometricCorrFromLatLon(lat, lon);
            end
            if nargin < 5
                time = GPS_Time();
            end
            
            N = GPS_SS.ELL_A ./ sqrt(1 - GPS_SS.ELL_E.^2 * sin(lat).^2);
            
            x = (N + h_ellips) .* cos(lon) .* cos(lat);
            y = (N + h_ellips) .* sin(lon) .* cos(lat);
            z = (N * (1 - GPS_SS.ELL_E.^2) + h_ellips) .* sin(lat);

            this = Coordinates.fromXYZ(x, y, z, time);
        end
    
        function this = fromCooFile(file_name)
            % Importing from a coo file XYZ and timestamp to a Coordinate
            % object
            this = Coordinates;
            if exist(file_name, 'file') == 2
                % Read and append
                [txt, lim] = Core_Utils.readTextFile(file_name, 3);
                if isempty(lim)
                    f = false;
                    timestamp = [];
                else
                    % Verify the file version (it should match 1.0):
                    id_ver = find(txt(lim(:,1) + 1) == 'F'); % +FileVersion
                    file_ok = not(isempty(regexp(txt(lim(id_ver, 1):lim(id_ver, 2)), '(?<=FileVersion[ ]*: )1.0', 'once')));
                    
                    % Point Name
                    timestamp = [];
                    if file_ok
                        id_line = find(txt(lim(:,1) + 1) == 'M'); % +MonitoringPoint
                        if isempty(id_line)
                            file_ok = false;
                        else
                            this.name = regexp(txt(lim(id_line, 1):lim(id_line, 2)), '(?<=MonitoringPoint[ ]*: ).*', 'match', 'once');
                        end
                    end
                    
                    % Data column (at the moment set here manually)
                    data_col = [2, 3, 4] + 1; % x, y, z
                    
                    % Data should be present
                    timestamp = [];
                    if file_ok
                        id_len_ok = find(lim(:,3)+1 >= 9);
                        data_start = id_len_ok(find(txt(lim(id_len_ok,1) + 9) == 't') + 1); % +DataStart
                        id_len_ok = find(lim(:,3)+1 >= 8);
                        data_stop = id_len_ok(find(txt(lim(id_len_ok,1) + 7) == 'd') -1); % +DataStop
                        if isempty(data_stop)
                            data_stop = size(lim, 1);
                        end
                        if isempty(data_start)
                            file_ok = false;
                        else
                            id_data = lim(data_start:data_stop,1);
                            % Read old timestamps
                            timestamp = datenum(txt(repmat(id_data, 1, 19) + repmat(0:18, numel(id_data), 1)), 'yyyy-mm-dd HH:MM:SS');
                        end
                    end
                    
                    % Import XYZ and time
                    if file_ok
                        this.xyz = nan(data_stop - data_start + 1, 3);
                        for l = 0 : (data_stop - data_start)
                            data_line = strsplit(txt(lim(data_start + l, 1) : lim(data_start + l, 2)), ';');
                            this.xyz(l + 1, 1) = str2double(data_line{data_col(1)});
                            this.xyz(l + 1, 2) = str2double(data_line{data_col(2)});
                            this.xyz(l + 1, 3) = str2double(data_line{data_col(3)});
                        end
                        this.time = GPS_Time(timestamp);
                    end
                    
                    % Check description field
                    % use the old one for the file
                end
            else
                log = Core.getLogger();
                log.addError(sprintf('%s cannot be imported', file_name));
            end
        end
    end
    
    % =========================================================================
    %    SHOWs
    % =========================================================================
    
    methods (Access = 'public')
        function fh = showPositionENU(coo_list)
            % Plot East North Up coordinates
            %
            % SYNTAX 
            %   this.showPositionENU(coo_list);
            fh = showCoordinatesENU(coo_list);
        end
        
        function fh_list = showCoordinates(mode, coo_list, coo_ref, n_obs)
            % Plot ENU or XYZ coordinates
            %
            % SYNTAX
            %   this.showCoordinates(coo_list);
            
            if strcmpi(mode, 'XYZ')
                mode = 'XYZ';
                axis_label = {'ECEF X', 'ECEF Y', 'ECEF Z'};
            else
                mode = 'ENU';
                axis_label = {'East', 'North', 'Up'};
            end
            fh_list = [];
            thr = 0.8;
            
            str_title{2} = sprintf('STD (detrended)');
            str_title{3} = sprintf('STD (detrended)');
            log = Core.getLogger();
            for i = 1 : numel(coo_list)
                pos = coo_list(i);
                
                if ~pos.isEmpty
                    if nargin < 3 || isempty(coo_ref)
                        if not(isempty(pos.name))
                            str_title{1} = sprintf('%s\nPosition stability %s [mm]\nSTD (detrended)', pos.name, mode);
                            fig_name = sprintf('d%s %s', mode, pos.name);
                        else
                            str_title{1} = sprintf('Position stability %s [mm]\nSTD (detrended)', mode);
                            fig_name = sprintf('d%s', mode);
                        end
                        
                        if strcmpi(mode, 'XYZ')
                            pos_diff = (pos.getXYZ - pos.getMedianPos.getXYZ) * 1e3;
                        elseif strcmpi(mode, 'ENU')
                            pos_diff = pos.getLocal(pos.getMedianPos) * 1e3;
                        end
                        flag_time = true;
                        if isa(pos.time, 'GPS_Time') && ~pos.time.isEmpty
                            t = pos.time.getMatlabTime;
                            if numel(t) < size(pos_diff,1)
                                log.addWarning(sprintf('Coordinates are corrupted, it seems that there are more coordinates than times\n plotting only the positions with time'))
                                pos_diff = pos_diff(1:numel(t),:);
                            elseif numel(t) > size(pos_diff,1)
                                log.addWarning(sprintf('Coordinates are corrupted, it seems that there are more times than coordinates\n plotting only the first positions'))
                                t = t(1:size(pos_diff,1),:);
                            end
                        else
                            flag_time = false;
                            t = (1 : size(pos_diff, 1))';
                        end
                    elseif nargin > 2 && not(isempty(coo_ref)) % plot baseline
                        if not(isempty(pos.name))
                            str_title{1} = sprintf('%s - %s\nBaseline stability %s [mm]\nSTD (detrended)', pos.name, coo_ref.name, mode);
                            fig_name = sprintf('d%s %s - %s', mode, pos.name, coo_ref.name);
                        else
                            str_title{1} = sprintf('Baseline stability %s [mm]\nSTD (detrended)', mode);
                            fig_name = sprintf('d%s bsl', mode);
                        end
                        
                        pos_diff = [];
                        if isa(pos.time, 'GPS_Time') && ~pos.time.isEmpty
                            [t_comm, idx_1, idx2] = intersect(round(coo_ref.time.getRefTime(pos.time.first.getMatlabTime)),round(pos.time.getRefTime(pos.time.first.getMatlabTime)));
                            t = pos.time.first.getMatlabTime + t_comm/86400;
                            if strcmpi(mode, 'XYZ')
                                pos_diff = (pos.xyz(idx2,:) - coo_ref.xyz(idx_1,:))*1e3;
                                pos_diff = bsxfun(@minus, pos_diff, median(pos_diff,1, 'omitnan'));
                            elseif strcmpi(mode, 'ENU')
                                pos_diff = Coordinates.cart2local(median(coo_ref.xyz,1,'omitnan'),pos.xyz(idx2,:) - coo_ref.xyz(idx_1,:) )*1e3;
                            end
                            flag_time = true;
                        else
                            if numel(coo_ref.xyz) == numel(pos.xyz)
                                pos_diff = Coordinates.cart2local(median(coo_ref.xyz,1,'omitnan'),pos.xyz - coo_ref.xyz)*1e3;
                                t = t(1:size(pos_diff,1),:);
                                flag_time = false;
                            else
                                log.addError(sprintf('No time in coordinates and number off coordinates in ref different from coordinate in the second receiver'))
                            end
                        end
                        pos_diff = bsxfun(@minus, pos_diff,median(pos_diff,1,'omitnan'));
                    end
                    
                    
                    fh = figure('Visible', 'off');  subplot(3,1,1); Core_UI.beautifyFig(fh);
                    fh.Name = sprintf('%03d: %s', fh.Number, fig_name);  fh.NumberTitle = 'off';
                    
                    if size(pos_diff, 1) > 1
                        if nargin >= 4 && n_obs > 0
                            id_ok = (max(1, size(pos_diff,1) - n_obs + 1)) : size(pos_diff,1);
                            pos_diff = pos_diff(id_ok, :);
                            t = t(id_ok);
                        end
                        
                        if numel(coo_list) == 1
                            color_order = Core_UI.getColor(1:3,3);
                        else
                            color_order = Core_UI.getColor(i * [1 1 1], numel(coo_list));
                        end      
                        
                        setAxis(fh, 1);
                        e = pos_diff(:,1);
                        if thr < 1
                            [data, lid_ko, trend] = strongFilterStaticData(e, 0.8, 7);
                            setAxis(fh, 1);
                            Core_Utils.plotSep(t, e, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); hold on;
                            Core_Utils.plotSep(t, data, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:));
                            e = data;
                        else
                            setAxis(fh, 1);
                            Core_Utils.plotSep(t, e, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(1,:)); hold on;
                            lid_ko = false(numel(t), 1);
                            trend = Core_Utils.interp1LS(t(~isnan(pos_diff(:,1))), pos_diff(~isnan(pos_diff(:,1)),1), 1, t);
                        end
                        setAxis(fh, 1);
                        ax(3) = gca(fh);
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        yl = minMax(e);
                        ylim([min(-20, yl(1)) max(20, yl(2))]);
                        if flag_time
                            setTimeTicks(4);
                        end
                        h = ylabel([axis_label{1} ' [mm]']); h.FontWeight = 'bold';
                        grid on;
                        str_title{1} = sprintf('%s %s%.2f', str_title{1}, iif(i>1, '- ', ''), std((e(~lid_ko) - trend(~lid_ko)), 'omitnan'));
                        h = title(str_title{1}, 'interpreter', 'none'); h.FontWeight = 'bold';
                        subplot(3,1,2);
                        
                        n = pos_diff(:,2);                        
                        if thr < 1
                            [data, lid_ko, trend] = strongFilterStaticData(n, 0.8, 7);
                            setAxis(fh, 2);
                         Core_Utils.plotSep(t, n, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); hold on;
                            Core_Utils.plotSep(t, data, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:)); 
                            n = data;
                        else
                            setAxis(fh, 2);
                        Core_Utils.plotSep(t, n, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(2,:)); hold on;
                            trend = Core_Utils.interp1LS(t(~isnan(pos_diff(:,2))), pos_diff(~isnan(pos_diff(:,2)),2), 1, t);
                        end
                        ax(2) = gca(fh);
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        yl = minMax(n);
                        ylim([min(-20, yl(1)) max(20, yl(2))]);
                        if flag_time
                            setTimeTicks(4);
                        end
                        h = ylabel([axis_label{2} ' [mm]']); h.FontWeight = 'bold';
                        str_title{2} = sprintf('%s %s%.2f', str_title{2}, iif(i>1, '- ', ''), std((n(~lid_ko) - trend(~lid_ko)), 'omitnan'));
                        h = title(str_title{2}, 'interpreter', 'none'); h.FontWeight = 'bold';
                        grid on;
                        subplot(3,1,3);
                        
                        up = pos_diff(:,3);                        
                        if thr < 1
                            [data, lid_ko, trend] = strongFilterStaticData(up, 0.8, 7);
                            setAxis(fh, 3);
                        Core_Utils.plotSep(t, up, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', [0.5 0.5 0.5]); hold on;
                            Core_Utils.plotSep(t, data, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:)); 
                            up = data;
                        else
                            trend = Core_Utils.interp1LS(t(~isnan(pos_diff(:,3))), pos_diff(~isnan(pos_diff(:,3)),3), 1, t);
                            setAxis(fh, 3);
                        Core_Utils.plotSep(t, up, '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:)); hold on;
                        end
                        ax(1) = gca(fh);
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        yl = minMax(up);
                        ylim([min(-20, yl(1)) max(20, yl(2))]);
                        if flag_time
                            setTimeTicks(4);
                        end
                        h = ylabel([axis_label{3} ' [mm]']); h.FontWeight = 'bold';
                        str_title{3} = sprintf('%s %s%.2f', str_title{3}, iif(i>1, '- ', ''), std((up(~lid_ko) - trend(~lid_ko)), 'omitnan'));
                        h = title(str_title{3}, 'interpreter', 'none'); h.FontWeight = 'bold';
                        grid on;
                        linkaxes(ax, 'x');
                        grid on;
                    else
                        log.addMessage('Plotting a single point static coordinates is not yet supported');
                    end
                    fh_list = [fh_list fh];
                    Core_UI.beautifyFig(fh);
                    Core_UI.addBeautifyMenu(fh);
                    fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); drawnow;
                end
            end
            
        end

        function [lid_ko, time, lid_ko_enu, trend_enu] = getBadSession(coo_list, coo_ref, n_obs)
            % Get outliers in East North Up coordinates
            %
            % SYNTAX 
            %   [lid_ko lid_ko_enu, trend_enu] = getBadSession(coo_list, coo_ref, n_obs);
            
            thr = 0.8;
            time = {};
            log = Core.getLogger();
            fh = figure('Visible', 'off'); Core_UI.beautifyFig(fh);
            for i = 1 : numel(coo_list)
                pos = coo_list(i);
                if ~pos.isEmpty
                    
                    if nargin == 1 || isempty(coo_ref)
                        enu_diff = pos.getLocal(pos.getMedianPos) * 1e3;
                        flag_time = true;
                        if isa(pos.time, 'GPS_Time') && ~pos.time.isEmpty
                            t = pos.time.getMatlabTime;
                            if numel(t) < size(enu_diff,1)
                                log.addWarning(sprintf('Coordinates are corrupted, it seems that there are more coordinates than times\n plotting only the positions with time'))
                                enu_diff = enu_diff(1:numel(t),:);
                            elseif numel(t) > size(enu_diff,1)
                                log.addWarning(sprintf('Coordinates are corrupted, it seems that there are more times than coordinates\n plotting only the first positions'))
                                t = t(1:size(enu_diff,1),:);
                            end
                        else
                            flag_time = false;
                            t = (1 : size(enu_diff, 1))';
                        end
                    elseif nargin > 1 % plot baseline
                        enu_diff = [];
                        if isa(pos.time, 'GPS_Time') && ~pos.time.isEmpty
                            [t_comm, idx_1, idx2] = intersect(round(coo_ref.time.getRefTime(pos.time.first.getMatlabTime)),round(pos.time.getRefTime(pos.time.first.getMatlabTime)));
                            t = pos.time.first.getMatlabTime + t_comm/86400;
                            enu_diff = Coordinates.cart2local(median(coo_ref.xyz,1,'omitnan'),pos.xyz(idx2,:) - coo_ref.xyz(idx_1,:) )*1e3;
                            flag_time = true;
                        else
                            if numel(coo_ref.xyz) == numel(pos.xyz)
                                enu_diff = Coordinates.cart2local(median(coo_ref.xyz,1,'omitnan'),pos.xyz - coo_ref.xyz)*1e3;
                                t = t(1:size(enu_diff,1),:);
                                flag_time = false;
                            else
                                log.addError(sprintf('No time in coordinates and number off coordinates in ref different from coordinate in the second receiver'))
                            end
                        end
                        enu_diff = bsxfun(@minus, enu_diff,median(enu_diff,1,'omitnan'));
                    end
                    
                    if nargin >= 2 && n_obs > 0
                        id_ok = (max(1, size(enu_diff,1) - n_obs + 1)) : size(enu_diff,1);
                        enu_diff = enu_diff(id_ok, :);
                        t = t(id_ok);
                    end
                    
                    e = enu_diff(1:numel(t),1);
                    [data, lid_ko_enu{i}(:,1), trend_enu{i}(:,1)] = strongFilterStaticData(e, 0.8, 7);
                    
                    n = enu_diff(:,2);
                    [data, lid_ko_enu{i}(:,2), trend_enu{i}(:,1)] = strongFilterStaticData(n, 0.8, 7);
                    
                    up = enu_diff(:,3);
                    [data, lid_ko_enu{i}(:,3), trend_enu{i}(:,1)] = strongFilterStaticData(up, 0.8, 7);
                    
                    lid_ko = isnan(e) | lid_ko_enu{i}(:,1) | isnan(n) | lid_ko_enu{i}(:,2) | isnan(up) | lid_ko_enu{i}(:,3);
                    time{i} = t;
                end
            end
            if numel(lid_ko) == 1
                lid_ko = lid_ko{1};
                lid_ko_enu = lid_ko_enu{1};
                trend_enu = trend_enu{1};
                time = time{1};
            end
        end

        function fh = showCoordinatesENU(coo_list, coo_ref, n_obs)
            % Plot East North Up coordinates
            %
            % SYNTAX 
            %   this.showCoordinatesENU(coo_list);
            
            switch nargin
                case 1, fh = showCoordinates('ENU', coo_list);
                case 2, fh = showCoordinates('ENU', coo_list, coo_ref);
                case 3, fh = showCoordinates('ENU', coo_list, coo_ref, n_obs);
            end
        end
        
        function fh = showPositionXYZ(coo_list)
            % Plot X Y Z coordinates
            %
            % SYNTAX
            %   this.showPositionXYZ(coo_list);
           fh = showCoordinatesXYZ(coo_list);
        end
        
        function fh = showCoordinatesXYZ(coo_list, coo_ref, n_obs)
            % Plot X Y Z coordinates
            %
            % SYNTAX
            %   this.showCoordinatesXYZ(coo_list);
            switch nargin
                case 1, fh = showCoordinates('XYZ', coo_list);
                case 2, fh = showCoordinates('XYZ', coo_list, coo_ref);
                case 3, fh = showCoordinates('XYZ', coo_list, coo_ref, n_obs);
            end
        end
        
        function fh = showPositionPlanarUp(coo_list)
            % Plot East North Up coordinates
            %
            % SYNTAX 
            %   this.showPositionENU(coo_list);
            fh = showCoordinatesENU(coo_list);
        end
        
        function fh = showCoordinatesPlanarUp(coo_list, coo_ref, n_obs)
            % Plot East North Up coordinates
            %
            % SYNTAX 
            %   this.showCoordinatesENU(coo_list);
            
            log = Core.getLogger();
            for i = 1 : numel(coo_list)
                pos = coo_list(i);
                if ~pos.isEmpty
                    
                    if nargin == 1 || isempty(coo_ref)
                        enu_diff = pos.getLocal(pos.getMedianPos) * 1e3;
                        flag_time = true;
                        if isa(pos.time, 'GPS_Time') && ~pos.time.isEmpty
                            t = pos.time.getMatlabTime;
                            if numel(t) < size(enu_diff,1)
                                log.addWarning(sprintf('Coordinates are corrupted, it seems that there are more coordinates than times\n plotting only the positions with time'))
                                enu_diff = enu_diff(1:numel(t),:);
                            elseif numel(t) > size(enu_diff,1)
                                log.addWarning(sprintf('Coordinates are corrupted, it seems that there are more times than coordinates\n plotting only the first positions'))
                                t = t(1:size(enu_diff,1),:);
                            end
                        else
                            flag_time = false;
                            t = (1 : size(enu_diff, 1))';
                        end
                    elseif nargin > 1 % plot baseline
                        enu_diff = [];
                        if isa(pos.time, 'GPS_Time') && ~pos.time.isEmpty
                            [t_comm, idx_1, idx2] = intersect(round(coo_ref.time.getRefTime(pos.time.first.getMatlabTime)),round(pos.time.getRefTime(pos.time.first.getMatlabTime)));
                            t = pos.time.first.getMatlabTime + t_comm/86400;
                            enu_diff = Coordinates.cart2local(median(coo_ref.xyz,1,'omitnan'),pos.xyz(idx2,:) - coo_ref.xyz(idx_1,:) )*1e3;
                            flag_time = true;
                        else
                            if numel(coo_ref.xyz) == numel(pos.xyz)
                                enu_diff = Coordinates.cart2local(median(coo_ref.xyz,1,'omitnan'),pos.xyz - coo_ref.xyz)*1e3;
                                t = t(1:size(enu_diff,1),:);
                                flag_time = false;
                            else
                                log.addError(sprintf('No time in coordinates and number off coordinates in ref different from coordinate in the second receiver'))
                            end
                        end
                        enu_diff = bsxfun(@minus, enu_diff,median(enu_diff,1,'omitnan'));
                    end
                    if size(enu_diff, 1) > 1
                        
                        if nargin >= 2 && n_obs > 0
                            id_ok = (max(1, size(enu_diff,1) - n_obs + 1)) : size(enu_diff,1);
                            enu_diff = enu_diff(id_ok, :);
                            t = t(id_ok);
                        end
                        str_title{1} = sprintf('Position detrended planar Up [mm]\nSTD (detrended)');
                        str_title{3} = sprintf('STD (detrended)');
                        fh = figure('Visible', 'off'); Core_UI.beautifyFig(fh);
                        if numel(coo_list) > 1
                            fh.Name = sprintf('%03d: dPUP MR', fh.Number); fh.NumberTitle = 'off';
                        else
                            fh.Name = sprintf('%03d: dPUP', fh.Number); fh.NumberTitle = 'off';
                        end
                        
                        if numel(coo_list) == 1
                            color_order = Core_UI.getColor(1:3,3);
                        else
                            color_order = Core_UI.getColor(i * [1 1 1], numel(coo_list));
                        end    
                        
                        trend_e = Core_Utils.interp1LS(t(~isnan(enu_diff(:,1))), enu_diff(~isnan(enu_diff(:,1)),1), 1, t);
                        trend_n = Core_Utils.interp1LS(t(~isnan(enu_diff(:,2))), enu_diff(~isnan(enu_diff(:,2)),2), 2, t);
                        trend_u = Core_Utils.interp1LS(t(~isnan(enu_diff(:,3))), enu_diff(~isnan(enu_diff(:,3)),3), 3, t);
                        enu_diff(:,1) = enu_diff(:,1) - trend_e(:);
                        enu_diff(:,2) = enu_diff(:,2) - trend_n(:);
                        enu_diff(:,3) = enu_diff(:,3) - trend_u(:);
                        
                        main_vb = uix.VBox('Parent', fh, ...
                            'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                        
                        tmp_box1 = uix.VBox('Parent', main_vb, ...
                            'Padding', 5, ...
                            'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                        tmp_box2 = uix.VBox('Parent', main_vb, ...
                            'Padding', 5, ...
                            'BackgroundColor', Core_UI.LIGHT_GREY_BG);
                        main_vb.Heights = [-2 -1];
                        Core_UI.beautifyFig(fh);
                        fh.Visible = 'on';
                        drawnow
                        fh.Visible = 'off';
                        ax = axes('Parent', tmp_box1);

                        % Plot parallel
                        max_e = ceil(max(abs(minMax(enu_diff)))/5) * 5;
                        max_n = ceil(max(abs(minMax(enu_diff)))/5) * 5;
                        max_r = ceil(sqrt(max_e^2 + max_n^2) / 5) * 5;
                       
                        % Plot circles of precision
                        az_l = 0 : pi/200: 2*pi;
                        % dashed
                        id_dashed = serialize(bsxfun(@plus, repmat((0:20:395)',1,5), (1:5)));
                        az_l(id_dashed) = nan;
                        step = max(1, floor((max_r/10)/5)) * 10;
                        if step > 10
                            step = round(step/50)*50;
                        end
                        if step > 100
                            step = round(step/100)*100;
                        end    
                        decl_s = ((step : step  : max_r));
                        for d = decl_s
                            x = cos(az_l).*d;
                            y = sin(az_l).*d;
                            plot(x,y,'color',[0.6 0.6 0.6], 'LineWidth', 2); hold on;
                            x = cos(az_l).*(d-step/2);
                            y = sin(az_l).*(d-step/2);
                            plot(x,y,'color',[0.75 0.75 0.75], 'LineWidth', 2); hold on;
                        end
                        
                        plot(enu_diff(:,1) + trend_e(:), enu_diff(:,2) + trend_n(:), 'o', 'MarkerSize', 4, 'LineWidth', 2, 'Color', color_order(1,:)); hold on;
                        axis equal;
                        h = ylabel('East [mm]'); h.FontWeight = 'bold';
                        h = xlabel('North [mm]'); h.FontWeight = 'bold';
                        ylim(max_r * [-1 1]);
                        xlim(max_r * [-1 1]);
                        grid on;
                        
                        str_title{1} = sprintf('%s %.2f - %.2f - %.2f', str_title{1}, std(enu_diff(:,1), 'omitnan'), std(enu_diff(:,2), 'omitnan'), std(enu_diff(:,3), 'omitnan'));
                        h = title(str_title{1}, 'interpreter', 'none'); h.FontWeight = 'bold';
                        h.FontWeight = 'bold';
                        
                        set(0, 'CurrentFigure', fh);;
                        ax = axes('Parent', tmp_box2);
                        up = enu_diff(:,3);                        
                        Core_Utils.plotSep(t, up + trend_u(:), '.-', 'MarkerSize', 15, 'LineWidth', 2, 'Color', color_order(3,:)); hold on;
                        ax(1) = gca(fh);
                        if (t(end) > t(1))
                            xlim([t(1) t(end)]);
                        end
                        yl = [-1 1] * max(abs([perc(up + trend_u(:), 0.05)*3 perc(up + trend_u(:), 0.95)*3]));
                        ylim([min(-20, yl(1)) max(20, yl(2))]);
                        setTimeTicks(4); h = ylabel('Up [mm]'); h.FontWeight = 'bold';
                        grid on;
                        drawnow;                        
                        Core_UI.beautifyFig(fh);
                        Core_UI.addBeautifyMenu(fh);
                        fh.Visible = iif(Core_UI.isHideFig, 'off', 'on'); 
                    else
                        log.addMessage('Plotting a single point static coordinates is not yet supported');
                    end
                end
            end            
        end
    end
    
    % =========================================================================
    %    OPERATIONS
    % =========================================================================
    
    methods (Access = 'public')
        function this = importCoo(this, file_name)
            this = Coordinates.fromCooFile(file_name);
        end
        
        function res = eq(coo1, coo2)
            %%% DESCRIPTION: check if two coordinates are equal
            d = sqrt(sum((coo1.xyz - coo2.xyz).^2, 2));
            res = d < coo1.precision;
        end
        
        function coo_diff = minus(coo1, coo2)
            coo_diff = struct('time', [], 'enu_diff', [], 'xyz_diff', []);
            if coo1.time.isEmpty || coo2.time.isEmpty
                coo_diff.time = GPS_Time;
                
                id_ok1 = 1 : min(size(coo1.xyz, 1), size(coo2.xyz, 1));
                id_ok2 = id_ok1;
            else
                [common_time, id_ok1, id_ok2] = intersect(coo1.time.getNominalTime.getMatlabTime, coo2.time.getNominalTime.getMatlabTime);
                coo_diff.time = GPS_Time(common_time);
            end
            
            coo_diff.enu_diff = coo1.getElement(id_ok1).getENU - coo2.getElement(id_ok2).getENU;
            coo_diff.xyz_diff = coo1.getElement(id_ok1).xyz - coo2.getElement(id_ok2).xyz;
            %coo_diff.enu_diff = Coordinates.cart2local(coo1.getElement(id_ok1).getMedianPos.xyz, coo_diff.xyz_diff);
        end
    end
        
    
    methods (Access = 'public', Static)
        function N = getOrthometricCorrFromLatLon(phi, lam, geoid, method)
            % SYNTAX:
            %   N = getOrthometricCorr(phi, lam, geoid);
            %
            % EXAMPLE:
            %   core = Core.getInstance;
            %   core.initGeoid();
            %   Coordinates.getOrthometricCorrFromLatLon(45.69 ./ 180*pi, 9.03 ./ 180*pi)
            %   % answer should be 46.1767008
            %
            % INPUT:
            %   phi     = geodetic latitude                [deg (rad only for legacy method)]
            %   lam     = geodetic longitude               [deg (rad only for legacy method)]
            %   geoid   = regular map in geocentric coordinates <default = EGM08 0.5x0.5 deg>
            %   method  = interpolation approach:
            %              - legacy
            %              - grid
            %              - grid_cubic  ( phi, lam are array of grid coordinates )
            %              - grid_akima  ( phi, lam are array of grid coordinates )
            %              - linear
            %              - natural
            %
            % OUTPUT:
            %   N       = geoid ondulation [m]
            %
            % DESCRIPTION:
            %   Get the geoid ondulation (orthometric correction)

            if (nargin < 3) || isempty(geoid)
                geoid = Core.getRefGeoid();
            end
            
            if (geoid.grid == 0)
                core = Core.getInstance(false);
                core.initGeoid();
                geoid = core.getRefGeoid();
            end
            
            if (nargin < 4) || isempty(method)
                method = 'legacy';
            end

            if  (geoid.ncols == 0 || geoid.nrows == 0)
                Core.initGeoid();
                geoid = Core.getRefGeoid();
            end
            N = zeros(numel(lam), 1);

            switch method
                case 'legacy'
                    for i = 1 : numel(lam)
                        N(i) = grid_bilin_interp(lam(i) / pi * 180, phi(i) / pi * 180, geoid.grid, geoid.ncols, geoid.nrows, geoid.cellsize, geoid.Xll, geoid.Yll, -9999);
                    end
                case 'grid'
                    x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
                    y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));

                    [xmg, ymg] = meshgrid(x_grid, y_grid);
                    N = interp2(xmg, ymg, geoid.grid, lam, phi, 'linear');
                case 'grid_cubic'
                    x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
                    y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));

                    [xmg, ymg] = meshgrid(x_grid, y_grid);
                    N = interp2(xmg, ymg, geoid.grid, lam, phi, 'cubic');
                case 'grid_akima'
                    x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
                    y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));

                    [xmg, ymg] = meshgrid(x_grid, y_grid);
                    N = interp2(xmg, ymg, geoid.grid, lam, phi, 'makima');
                case 'linear'
                    x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
                    y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));
                    [xmg, ymg] = meshgrid(x_grid, y_grid);
                    finterp = scatteredInterpolant(xmg(:), ymg(:), geoid.grid(:), 'linear');

                    N = finterp(lam, phi);
                case 'natural'
                    x_grid = geoid.Xll + geoid.cellsize * (0 : geoid.ncols - 1);
                    y_grid = fliplr(geoid.Yll + geoid.cellsize * (0 : geoid.nrows - 1));

                    [xmg, ymg] = meshgrid(x_grid, y_grid);
                    finterp = scatteredInterpolant(xmg(:), ymg(:), geoid.grid(:), 'natural');
                    N = finterp(lam, phi);
            end
        end
        
        function [lat, lon, h, lat_geoc] = cart2geod(xyz, y, z)
            % Get Geodetic coordinates from xyz cartesian coordinates
            %
            % SYNTAX:
            %   [lat, lon, h, lat_geoc] = Coordinates.cart2geod(xyz)
            %   [lat, lon, h, lat_geoc] = Coordinates.cart2geod(x, y, z)
            
            if nargin == 1
                z = xyz(:, 3);
                y = xyz(:, 2);
                xyz = xyz(:, 1);
            end
            
            a = Coordinates.ELL_A;
            e = Coordinates.ELL_E;
            
            % radius computation
            r = sqrt(xyz.^2 + y.^2 + z.^2);
            
            % longitude
            lon = atan2(y,xyz);
            
            % geocentric latitude
            lat_geoc = atan(z./sqrt(xyz.^2 + y.^2));
            
            % Coordinates transformation
            psi = atan(tan(lat_geoc)/sqrt(1-e^2));
            
            lat = atan((r.*sin(lat_geoc) + e^2*a/sqrt(1-e^2) * (sin(psi)).^3) ./ ...
                (r.*cos(lat_geoc) - e^2*a * (cos(psi)).^3));
            
            if nargout > 2
                N = a ./ sqrt(1 - e^2 * sin(lat).^2);
                
                % height
                h = r .* cos(lat_geoc)./cos(lat) - N;
            end
        end
        
        function [lat_geoc, lon, r] = cart2geoc(xyz, y, z)
            % Get Geocentric sphericals coordinates from xyz cartesian coordinates
            %
            % SYNTAX:
            %   [lat_geoc, lon, r] = Coordinates.cart2geoc(xyz)
            %   [lat_geoc, lon, r] = Coordinates.cart2geoc(x, y, z)
            
            if nargin == 1
                z = xyz(:, 3);
                y = xyz(:, 2);
                xyz = xyz(:, 1);
            end
            
            a = Coordinates.ELL_A;
            e = Coordinates.ELL_E;
            
            % radius computation
            r = sqrt(xyz.^2 + y.^2 + z.^2);
            
            % longitude
            lon = atan2(y,xyz);
            
            % geocentric latitude
            lat_geoc = atan(z./sqrt(xyz.^2 + y.^2));            
        end
        
        function [east, north, utm_zone] = geod2plan(lat, lon)
            % Conversion from geodetic coordinates to planimetric coordinates (UTM WGS84).
            %
            % SYNTAX:
            %   [east, north, utm_zone] = geod2plan(lat, lon);
            %
            % INPUT:
            %   lat = latitude [rad]
            %   lon = longitude [rad]
            %
            % OUTPUT:
            %   east     = east Coordinates [m]
            %   north    = north Coordinates [m]
            %   utm_zone = UTM zone
            %       
            
            %number of input points
            n = size(lat);
            n = n(1,1);
            
            if (n > 0)
                %pre-allocation
                M = zeros(n,1);
                north_south = zeros(n,1);
                utm_zone(n,:) = '60 X';
                north = zeros(n,1);
                east = zeros(n,1);
                
                % conversion algorithm
                % UTM parameters
                EsMc = 500000; % false East
                
                % WGS84 ellipsoid parameters
                % semi-major equatorial axis [m]
                ell_a = Coordinates.ELL_A;
                
                % squared eccentricity (a^2-b^2)/a^2
                e2 = Coordinates.E2;
                
                % contraction factor
                contr = 0.9996;
                
                ell_a = ell_a * contr;
                ecc4 = e2 * e2;
                ecc6 = ecc4 * e2;
                ecc8 = ecc6 * e2;
                
                k0 = ell_a * (e2 / 4 + ecc4 * 3 / 64 + ecc6 * 5 / 256 + ecc8 * 175 / 16384);
                k = ell_a - k0;
                k1 = ell_a * (ecc4 * 13 / 96 + ecc6 * 59 / 384 + ecc8 * 1307 / 8192);
                k2 = ell_a * (ecc6 * 61 / 484 + ecc8 * 609 / 2048);
                k3 = ell_a * (ecc8 * 49561 / 322560);
                c1 = (e2 * 5 - ecc4) / 6;
                c2 = (ecc4 * 104 - ecc6 * 45) / 120;
                c3 = ecc6 * 1237 / 1260;
                
                % Sines, cosines and latitude powers
                Elix = [lat lon];
                latsessadec = lat ./ pi .* 180;
                lonsessadec = lon ./ pi .* 180;
                
                fiSin(:,1) = sin(Elix(:,1));
                fiCos(:,1) = cos(Elix(:,1));
                fiSin2(:,1) = fiSin .* fiSin;
                fiSin4(:,1) = fiSin2 .* fiSin2;
                fiSin6(:,1) = fiSin4 .* fiSin2;
                
                % UTM zone finding
                for i = 1 : n
                    
                    M(i, 1) = fix((180 + lonsessadec(i, 1)) / 6) + 1;
                    
                    if latsessadec(i, 1) >= 0
                        north_south(i, 1) = 1; %1 north, 0 south
                    else
                        north_south(i, 1) = 0;
                    end
                    
                    if     (latsessadec(i, 1) < -72), letter = 'C';
                    elseif (latsessadec(i, 1) < -64), letter = 'D';
                    elseif (latsessadec(i, 1) < -56), letter = 'E';
                    elseif (latsessadec(i, 1) < -48), letter = 'F';
                    elseif (latsessadec(i, 1) < -40), letter = 'G';
                    elseif (latsessadec(i, 1) < -32), letter = 'H';
                    elseif (latsessadec(i, 1) < -24), letter = 'J';
                    elseif (latsessadec(i, 1) < -16), letter = 'K';
                    elseif (latsessadec(i, 1) <  -8), letter = 'L';
                    elseif (latsessadec(i, 1) <   0), letter = 'M';
                    elseif (latsessadec(i, 1) <   8), letter = 'N';
                    elseif (latsessadec(i, 1) <  16), letter = 'P';
                    elseif (latsessadec(i, 1) <  24), letter = 'Q';
                    elseif (latsessadec(i, 1) <  32), letter = 'R';
                    elseif (latsessadec(i, 1) <  40), letter = 'S';
                    elseif (latsessadec(i, 1) <  48), letter = 'T';
                    elseif (latsessadec(i, 1) <  56), letter = 'U';
                    elseif (latsessadec(i, 1) <  64), letter = 'V';
                    elseif (latsessadec(i, 1) <  72), letter = 'W';
                    else,                             letter = 'X';
                    end
                    
                    utm_zone(i,:) = sprintf('%02d %c', nan2zero(M(i, 1)), letter);
                end
                
                for i = 1 : n
                    % Longitude of the central meridian
                    LonMeridianoCentrale = -177 + 6 * (M(i, 1) - 1);
                    
                    % Distance of the point from the central meridian
                    % la_sd --> distance in decimal degrees
                    % la    --> distance in radians
                    la_sd = lonsessadec(i, 1) - LonMeridianoCentrale;
                    la = la_sd / 180 * pi;
                    
                    if la == 0
                        laSin = 0;
                        laCos = 1;
                        laCot = 0;
                    else
                        laSin = sin(la);
                        laCos = cos(la);
                        laCot = laCos / laSin ;
                    end
                    
                    % longitude with respect to central meridian
                    laCot2 = laCot * laCot;
                    
                    % psi
                    psi = Elix(i, 1) - e2 * fiSin(i, 1) * fiCos(i, 1) *(1 + c1 * fiSin2(i, 1) + c2 * fiSin4(i, 1) + c3 * fiSin6(i, 1));
                    psiSin = sin(psi);
                    psiCos = cos(psi);
                    psiTan = psiSin / psiCos ;
                    psiSin2 = psiSin * psiSin ;
                    
                    % omega
                    ome = atan(psiTan / laCos);
                    
                    % sigma
                    if laSin ~= 0
                        sigSin =laSin * psiCos;
                        sig = asin(sigSin);
                    else
                        sigSin = 0;
                        sig = 0;
                    end
                    
                    sigSin2 = sigSin * sigSin;
                    sigSin4 = sigSin2 * sigSin2;
                    sigCos2 = 1 - sigSin2;
                    sigCos4 = sigCos2 * sigCos2;
                    
                    % chi
                    chi = sig / 2 + pi / 4;
                    chiTan = tan(chi);
                    chiLog = log(chiTan);
                    
                    % constants
                    aa = psiSin * psiCos * laCos * (1 + sigSin2) / sigCos4;
                    bb = sigSin * (sigCos2 - 2 * psiSin2) / sigCos4;
                    
                    if laCot ~= 0
                        a1 = (psiSin2 - sigSin4 * laCot2)/ sigCos4;
                        b1 = 2 * sigSin2 * psiSin * laCot / sigCos4;
                    else
                        a1 = psiSin2 / sigCos4 ;
                        b1 = 0;
                    end
                    a2 = a1 * a1 - b1 * b1 ;
                    b2 = 2 * a1 * b1 ;
                    a3 = a1 * a2 - b1 * b2 ;
                    b3 = a1 * b2 + b1 * a2 ;
                    rr = k0 - a1 * k1 + a2 * k2 - a3 * k3;
                    tt = b1 * k1 - b2 * k2 + b3 * k3;
                    
                    % X/Y coordinates
                    xx = k * ome + aa * rr + bb * tt;
                    yy = k * chiLog + bb * rr - aa * tt;
                    
                    % North and East
                    if north_south(i, 1) == 1
                        NoEq = 0;
                    else
                        NoEq = 10000000;
                    end
                    
                    north(i, 1) = NoEq + xx;
                    east(i, 1) = EsMc + yy;
                end                
            else
                utm_zone = [];
                north = [];
                east = [];
            end
        end
        
        function [loc,  rot_mat] = cart2local(xyz_ref, xyz_baseline)
            % cart2local: from geocentric cartesian baselines (DX) to local coordinates in X0
            %
            % SYNTAX
            %   [loc,  rot_mat] = Coordinates.cart2local(xyz_ref, baseline)
            [lat, lon] = Coordinates.cart2geod(xyz_ref);
            rot_mat = [ -sin(lon) cos(lon) 0; -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)];
            loc = (rot_mat * xyz_baseline')';
        end
        
        function [xyz, rot_mat] = local2cart(xyz_ref, local_baseline)
            % local2cart: from local coordinates in X0 to geocentric cartesian baselines (DX)
            %
            % SYNTAX
            %   [loc,  rot_mat] = Coordinates.local2cart(xyz_ref, local_baseline)
            [lat, lon] = Coordinates.cart2geod(xyz_ref);
            rot_mat = [ -sin(lon) cos(lon) 0; -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)]';
            xyz = (rot_mat * local_baseline')';
        end
        
        function [loc,  rot_mat] = cart2loca(xyz_ref, xyz_baseline)
            [loc,  rot_mat] = Coordinates.cart2local(xyz_ref, xyz_baseline); % maybe the function was used with the wrong name
        end
        
        function [xyz_baseline, rot_mat] = loca2cart(xyz_ref, loc)
            % loca2cart: from local coordinates in X0 to geocentric cartesian baselines (DX)
            %
            % SYNTAX
            %   [baseline, rot_mat] = Coordinates.loca2cart(xyz_ref, loca)
            [lat, lon] = Coordinates.cart2geod(xyz_ref);
            rot_mat = [ -sin(lon) cos(lon) 0; -sin(lat)*cos(lon) -sin(lat)*sin(lon) cos(lat); cos(lat)*cos(lon) cos(lat)*sin(lon) sin(lat)]';
            xyz_baseline = (rot_mat * loc)';
        end        
                
        function fh = showCompareENU(coo_list)
            % Plot East North Up coordinates
            %
            % SYNTAX 
            %   fh = Coordinates.showCoordinatesENU(coo_list);
            
            fh = coo_list.showPositionENU();
        end
        
        function bslCompare(coo_list_0, coo_list_1, id_ref, spline_base, outlier_thr, y_lim)
            % Compare two sets of receivers coordinates (e.g. coming from two different execution)
            %
            % INPUT
            %   coo_list_0           first set of GNSS receivers or coordinates
            %   coo_list_1           second set of GNSS receivers or coordinates
            %   id_ref               id of the reference receiver
            %                        DEFAULT: last receiver
            %   spline_base          spline base in days (for removing splines)
            %                        DEFAULT: 365/4
            %   outlier_threshold    threshold on the baseline difference for the
            %                        spline and STD computation [mm]
            %                        DEFAULT: 2 mm
            %   y_lim                3 x 2 ylimits
            %
            % SINTAX
            %   Coordinates.bslCompare(rec_list_0, rec_list_1, id_ref, spline_base, outlier_thr)
            %
            % EXAMPLE
            %   Coordinates.bslCompare(nomp.core.rec, core.rec, 7, 365/4, 2);
            
            if nargin < 3 || isempty(id_ref)
                id_ref = numel(coo_list_0); % set to the last;
            end
            if nargin < 4 || isempty(spline_base)
                spline_base = 365/4;
            end
            % Set outlier treshold (on baseline differrence);
            if nargin < 5 || isempty(outlier_thr)
                outlier_thr = 2;
            end
            
            if isa(coo_list_0, 'Coordinates')
                coo_ref0 = coo_list_0(id_ref);
            else
                coo_ref0 = coo_list_0(id_ref).getPos;
            end
            if isa(coo_list_1, 'Coordinates')
                coo_ref1 = coo_list_1(id_ref);
            else
                coo_ref1 = coo_list_1(id_ref).getPos;
            end
            if nargin < 6 || isempty(y_lim)
                y_lim = [-5 5; -10 10; -10 10];
            elseif spline_base == 0
                spline_base = -1;
            end
            
            Core.getLogger.addMonoMessage('\nProcessing gain - solution 0 vs 1\n----------------------------------------');

            for r = setdiff(1 : numel(coo_list_1), id_ref)
                if isa(coo_list_0, 'Coordinates')
                    coo0 = coo_list_0(r);
                else
                    coo0 = coo_list_0(r).getPos;
                end
                if isa(coo_list_1, 'Coordinates')
                    coo1 = coo_list_1(r);
                else
                    coo1 = coo_list_1(r).getPos;
                end
                
                % Extract the baselines
                [t_comm, idx1, idx2] = intersect(round(coo_ref0.time.getRefTime(coo0.time.first.getMatlabTime)),round(coo0.time.getRefTime(coo0.time.first.getMatlabTime)));
                t0 = coo0.time.first.getMatlabTime + t_comm/86400;
                enu_diff0 = Coordinates.cart2local(median(coo_ref0.xyz,1,'omitnan'),coo0.xyz(idx2,:) - coo_ref0.xyz(idx1,:) )*1e3;
                id_ok = not(any(isnan(enu_diff0), 2)); % discard NaN
                t0 = t0(id_ok);
                enu_diff0 = enu_diff0(id_ok,:);
                
                [t_comm, idx1, idx2] = intersect(round(coo_ref1.time.getRefTime(coo1.time.first.getMatlabTime)),round(coo1.time.getRefTime(coo1.time.first.getMatlabTime)));
                t1 = coo1.time.first.getMatlabTime + t_comm/86400;
                enu_diff1 = Coordinates.cart2local(median(coo_ref1.xyz,1,'omitnan'),coo1.xyz(idx2,:) - coo_ref1.xyz(idx1,:) )*1e3;
                id_ok = not(any(isnan(enu_diff1), 2)); % discard NaN
                t1 = t1(id_ok);
                enu_diff1 = enu_diff1(id_ok,:);
                
                % Sync the baselines
                [t_comm, idx0, idx1] = intersect(round(t0*86400)/86400,round(t1*86400)/86400, 'stable');
                enu_diff = enu_diff0(idx0,:) - enu_diff1(idx1,:);
                id_ok = abs(enu_diff - median(enu_diff, 'omitnan')) < outlier_thr;
                
                fh = figure; Core_UI.beautifyFig(fh, 'dark'); drawnow
                % Plot the baseline difference ----------------------------------------
                subplot(3,1,1);
                tmp = bsxfun(@minus, enu_diff, [-20 0 20]);
                plotSep(t_comm, tmp, '.-', 'MarkerSize', 15, 'LineWidth', 2);
                for c = 1 : 3
                    std_enu(:,c) = std(enu_diff(id_ok(:,c), c), 'omitnan');
                end
                legend(sprintf('East (%.2f mm)', std_enu(1)), sprintf('North (%.2f mm)', std_enu(2)), sprintf('Up (%.2f mm)', std_enu(3)), 'location', 'EastOutside');
                ylim(y_lim(1,:));
                xlim(minMax(t_comm));
                setTimeTicks();
                ax(1) = gca;
                ylabel(sprintf('Baseline\ndifference'));
                grid minor
                if isempty(coo1.name)
                    trg_rec_name = sprintf('%d', r);
                else
                    trg_rec_name = coo1.name;
                end
                if isempty(coo_ref1.name)
                    ref_rec_name = sprintf('%d', id_ref);
                else
                    ref_rec_name = coo_ref1.name;
                end
                title(sprintf('Baseline %s - %s\\fontsize{5} \n', trg_rec_name, ref_rec_name), 'FontSize', 16);
                
                
                % Plot the baseline (filtered by spline) of the solution with no MP ---
                subplot(3,1,2);
                tmp = enu_diff0(idx0,:) - median(enu_diff0(idx0,:), 'omitnan');
                tmp0 = enu_diff0 - median(enu_diff0(idx0,:), 'omitnan');
                splined = zeros(size(tmp));
                splined0 = zeros(size(tmp0));
                if spline_base > 0
                    for c = 1 : 3
                        ttmp = t_comm(id_ok(:,c));
                        [filtered, ~, ~, splined(:,c)] = splinerMat(ttmp, tmp(id_ok(:,c),c), spline_base, 1e-8, t_comm);
                        [~, ~, ~, splined0(:,c)] = splinerMat(ttmp, tmp(id_ok(:,c),c), spline_base, 1e-8, t0);
                    end
                end
                tmp = tmp - splined; % remove splines
                tmp0 = tmp0 - splined0; % remove splines
                tmp0 = bsxfun(@minus, tmp0, [-20 0 20]);
                plotSep(t0, tmp0, '.-', 'MarkerSize', 15, 'LineWidth', 2);
                
                for c = 1 : 3
                    std_enu0(:,c) = std(tmp(id_ok(:,c), c), 'omitnan');
                end
                %plotSep(t_comm, splined, '.-', 'MarkerSize', 1, 'LineWidth', 1, 'Color', 'k');
                legend(sprintf('East (%.2f mm)', std_enu0(1)), sprintf('North (%.2f mm)', std_enu0(2)), sprintf('Up (%.2f mm)', std_enu0(3)), 'location', 'EastOutside');
                if spline_base ~= 0
                    ylim(y_lim(2,:));
                end
                xlim(minMax(t_comm));
                setTimeTicks();
                ax(2) = gca;
                ylabel(sprintf('Baseline 0\n(reduced)'));
                grid minor
                
                % Plot the baseline (filtered by spline) of the solution with MP ------
                subplot(3,1,3);
                tmp = enu_diff1(idx1,:) - median(enu_diff1(idx1,:), 'omitnan');
                tmp1 = enu_diff1 - median(enu_diff1, 'omitnan');
                splined = zeros(size(tmp));
                splined1 = zeros(size(tmp1));
                if spline_base > 0
                    for c = 1 : 3
                        ttmp = t_comm(id_ok(:,c));
                        [filtered, ~, ~, splined(:,c)] = splinerMat(ttmp, tmp(id_ok(:,c),c), spline_base, 1e-8, t_comm);
                        [~, ~, ~, splined1(:,c)] = splinerMat(ttmp, tmp(id_ok(:,c),c), spline_base, 1e-8, t1);
                    end
                end
                tmp = tmp - splined; % remove splines
                tmp1 = tmp1 - splined1; % remove splines
                tmp1 = bsxfun(@minus, tmp1, [-20 0 20]);
                plotSep(t1, tmp1, '.-', 'MarkerSize', 15, 'LineWidth', 2);
                
                for c = 1 : 3
                    std_enu1(:,c) = std(tmp(id_ok(:,c), c), 'omitnan');
                end
                legend(sprintf('East (%.2f mm)', std_enu1(1)), sprintf('North (%.2f mm)', std_enu1(2)), sprintf('Up (%.2f mm)', std_enu1(3)), 'location', 'EastOutside');
                if spline_base ~= 0
                    ylim(y_lim(3,:));
                end
                xlim(minMax(t_comm));
                setTimeTicks();
                ax(3) = gca;
                ylabel(sprintf('Baseline 1\n(reduced)'));
                grid minor
                
                Core_UI.beautifyFig(fh, 'dark');
                Core_UI.addExportMenu(fh);
                linkaxes(ax, 'x');
                
                Core.getLogger.addMonoMessage(sprintf('Baseline %d - %d) %5.2f %% %5.2f %% %5.2f %%', r, id_ref, (100*((std_enu0 - std_enu1) ./ std_enu0))));
            end
        end
    end
    
    % =========================================================================
    %    TESTS
    % =========================================================================
    
    methods (Static, Access = 'public')
        function test()
            % Testing function, tests some basic transformations
            % To be done
            %
            % SYNTAX
            %   test()
            log = Core.getLogger();
            
            log.addMessage('Testing Class Coordinates');
            tic;
            pos_diff = 0;

            if pos0 == pos1
                log.addStatusOk('Passed');
            else
                log.addWarning(sprintf('Difference greater than 0.2 ms: %e',t_diff));
            end
            toc
        end
    end
end
