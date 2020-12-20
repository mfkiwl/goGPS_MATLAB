function [data, lid_ko, trend] = strongFilterStaticData(data, robustness_perc, n_sigma)
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

    if nargin < 2
        robustness_perc = 0.8;
    end
    if nargin < 3
        n_sigma = 6;
    end
    [tmp, trend] = strongDeTrend(data, robustness_perc, 1-((1-robustness_perc)/2), n_sigma);
    thr = n_sigma * strongStd(tmp, robustness_perc);

    lid_ko = abs(tmp) > thr;
    % figure; plot(data, 'Color', [0.5 0.5 0.5]);
    data(lid_ko) = nan;
    % hold on; plot(data, '.-b', 'LineWidth', 2)
    % plot(tmp,'g');
end