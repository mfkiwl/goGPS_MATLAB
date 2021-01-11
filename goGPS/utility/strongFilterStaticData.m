function [data, lid_ko, trend, spline] = strongFilterStaticData(data, robustness_perc, n_sigma, spline_base)
% Returns the data removing outliers (spikes)
%
% INPUT:
%   data                column array of values
%   robustness_perc     maximum percentage of date with no outliers
%
% SYNTAX:
%   [data, id_ko] = strongFilterStaticData(data, robustness_perc)
    
%--- * --. --- --. .--. ... * ---------------------------------------------
%               ___ ___ ___
%     __ _ ___ / __| _ | __|
%    / _` / _ \ (_ |  _|__ \
%    \__, \___/\___|_| |___/
%    |___/                    v 1.0b8
%
%--------------------------------------------------------------------------
%  Copyright (C) 2020 Gatti Andrea, Giulio Tagliaferro, Eugenio Realini
%  Written by:       Andrea Gatti
%  Contributors:     Andrea Gatti ...
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

if any(data(:,end))
    if nargin < 2
        robustness_perc = 0.8;
    end
    if nargin < 3
        n_sigma = 6;
    end
    if nargin < 4 || numel(spline_base) ~= 2
        spline_base = [28, 2.5];
    end
    flag_time = false;
    idf = [];
    if size(data,2) >= 2
        if size(data,2) >= 3
            data_var = data(:,3);
        end
        time = data(:,1);
        rate = round(median(diff(time*86400)))/86400;
        time_full = linspace(time(1), time(end), round((time(end) - time(1)) / rate + 1))';
        [~, idf, idr] = intersect(round((time_full-rate/2)/rate), round((time-rate/2)/rate));
        tmp = data(:,2);
        data = nan(numel(time_full), 1);
        if numel(idr) < numel(tmp)
            Core.getLogger.addWarning('StrongFilter is loosing some observations out of sync');
        end
        data(idf) = tmp(idr);
        flag_time = 1;
    end
    [tmp, trend] = strongDeTrend(data, robustness_perc, 1-((1-robustness_perc)/2), n_sigma);
    
    if any(tmp) && flag_time && (numel(data(idf)) > 4)
        if (numel(tmp(idf)) > 11)
            if size(data,2) >= 3
                spline = splinerMat(time, [data(idf) data_var], spline_base(1), 1e-6); % one month splines
            else
                spline = splinerMat(time, [data(idf) tmp(idf).^2], spline_base(1), 1e-6); % one month splines
            end
            tmp(idf) = tmp(idf) - spline + trend(idf);
            spline = splinerMat(time, [data(idf) abs(tmp(idf))], spline_base(2), 1e-6); % one week splines
            spline = splinerMat(time, [data(idf) abs(data(idf) - spline)], spline_base(2), 1e-6); % one week splines
            spline = splinerMat(time, [data(idf) (data(idf) - spline).^2], spline_base(2), 1e-6); % one week splines
            tmp = data(idf) - spline;
        else
            tmp = tmp(idf);
            spline = trend(idf);
        end
    else
        if flag_time
            spline = trend(idf);
        else
            spline = trend;
        end
    end
    if flag_time
        data = data(idf);
        trend = trend(idf);
    end
    thr = n_sigma * strongStd(tmp, robustness_perc);
    lid_ko = abs(tmp) > thr;
    % figure; plot(data, 'Color', [0.5 0.5 0.5]);
    data(lid_ko) = nan;
    % hold on; plot(data, '.-b', 'LineWidth', 2)
    % plot(tmp,'g');
else
    % no data
    data = data(:,end);
    lid_ko = true(size(data));
    trend = data;
    spline = data;
end
end