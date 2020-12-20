function resetTimeTicks(h, num, format)
% reset function for setTimeTicks
% SEE ALSO: setTimeTicks

%--------------------------------------------------------------------------
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
%
%--------------------------------------------------------------------------
% 01100111 01101111 01000111 01010000 01010011
%--------------------------------------------------------------------------

    if (nargin == 1) 
        num = 4;
        format = 'yyyy/mm/dd HH:MM';
    end
    if isa(h, 'matlab.ui.Figure')
        h = get(h,'Children');
    end
    for i = 1 : numel(h)
        try
            h(i).TickLength = [0.005 0.01];
            ax = double(axis(h(i)));            

            step = (ax(2)-ax(1))/(num);
            time_span = ax(2)-ax(1);

            % depending on the time_span round ticks
            round_val = 1 / 86400; % one second
                        
            if time_span > 8                 % greater than 20days
                round_val =  1;               % round on 1 day
            elseif time_span > 4              % greater than 5 days
                round_val =  12 / 24;         % round on 12 hour
            elseif time_span > 1              % greater than a day
                round_val =  1 / 24;          % round on 1 hour
            elseif time_span > 12 / 24        % greater than 12 hours
                round_val =  30 / (24 * 60);  % round on 30 minutes
            elseif time_span > 5 / 24         % greater than 5 hours
                round_val =  15 / (24 * 60);  % round on 15 minutes
            elseif time_span > 1 / 24         % greater than 1 hour
                round_val =  5 / (24 * 60);   % round on 5 minutes
            elseif time_span > 10 / (24 * 60) % greater than 10 minutes
                round_val =  1 / (24 * 60);   % round on 1 minute
            elseif time_span > 5 / (24 * 60)  % greater than 5 minutes
                round_val =  30 / 86400;      % round on 30 seconds
            elseif time_span > 1 / (24 * 60)  % greater than 1 minutes
                round_val =  5 / 86400;       % round on 5 seconds
            end

            tick_pos = (ax(1) + (step/2)) : step : (ax(2) - (step/2));
            
            c_out = 0;
            if round_val >= 1 / (24 * 60) % rounding is greater than 1 minute
                c_out = 3;                % do not show seconds
            elseif round_val >= 1 / 24    % rounding is greater than 1 hour
                c_out = 6;                % do not show minutes
            elseif round_val >= 1         % greater than 1 day
                c_out = 9;                % do not show hours
            end
                
            tick_pos = round(tick_pos ./ round_val) .* round_val; % round ticks
            tick_pos((tick_pos < ax(1)) | (tick_pos > ax(2))) = []; % delete ticks outside figure;
            tick_pos = unique(tick_pos);
            
            if all(mod(tick_pos, 1) == 0)
                round_val = 1;
            end
            set(h(i), 'XTick', tick_pos);
            if strcmp(format, 'auto')
                last_date = [0 0 0 0 0 0];
                h.XTickLabel = cell(size(tick_pos));

                x_labels = cell(size(tick_pos));
                for l = 1 : numel(tick_pos)
                    str_time_format = 'HH:MM:SS';
                    new_date = datevec(double(tick_pos(l)));
                    if last_date(1) ~= new_date(1)        % if the year is different
                        str_time_format = [char(32 * ones(1, floor(c_out/2) + 5)) str_time_format(1:end - c_out) '-newlinemmm dd, yyyy'];
                    elseif last_date(2) ~= new_date(2)    % if the month is different 
                        str_time_format = [char(32 * ones(1, floor(c_out/2) + 5)) str_time_format(1:end - c_out) '-newlinemmm dd, yyyy'];
                    elseif last_date(3) ~= new_date(3)    % if the day is different 
                        str_time_format = [char(32 * ones(1, floor(c_out/2) + 5)) str_time_format(1:end - c_out) '-newlinemmm dd, yyyy'];
                    elseif last_date(4) ~= new_date(4)    % if the hour is different 
                        str_time_format = str_time_format(1:end - c_out);
                    else % if last_date(5) ~= new_date(5)    % if the hour is different
                        str_time_format = str_time_format(1:end - c_out);
                    end                    
                    if isempty(str_time_format)
                        h.XTickLabel{l} = '';
                    else
                        if round_val == 1
                            tmp = datestr(tick_pos(l), regexp(str_time_format, '(?<=line).*','match', 'once'));
                        else
                            tmp = datestr(tick_pos(l), str_time_format);
                            tmp(tmp == '-') = '\';
                        end
                        h.XTickLabel{l} = tmp;
                    end
                    last_date = new_date;
                end
            else
                datetick(h(i),'x',format,'keepticks');
            end
            axis(h(i),ax);
        catch ex
            Core_Utils.printEx(ex);
        end
    end
end
