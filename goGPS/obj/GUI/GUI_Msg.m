%   CLASS GUI_Msg
% =========================================================================
%
% DESCRIPTION
%   class to manage the Main Message Window of goGPS
%
% EXAMPLE
%   ui = GUI_Msg.getInstance();
%
% FOR A LIST OF CONSTANTs and METHODS use doc Core_UI


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

classdef GUI_Msg < GUI_Unique_Win   
    properties (Constant)
        WIN_NAME = 'goGPS_Msg_Win';
    end
    
    properties (Constant, Access = 'protected')
        BG_COLOR = Core_UI.DARK_GREY_BG;
    end
      
    %% PROPERTIES GUI
    % ==================================================================================================================================================
    properties
        win            % Handle to this window
        jedt           % j edit handle (java logger element)
    end    
        
    %% METHOD CREATOR
    % ==================================================================================================================================================
    methods  (Static, Access = private)
        function this = GUI_Msg()
            % GUI_MAIN object creator
            this.init();
            this.openGUI();
        end
    end         
    
    methods (Static, Access = public)
        function this = getInstance(flag_skip_init)
            % Get the persistent instance of the class
            
            persistent unique_instance_gui_msg__
            
            if isempty(unique_instance_gui_msg__) || ~ishandle(unique_instance_gui_msg__.win)
                this = GUI_Msg();
                unique_instance_gui_msg__ = this;
            else
                this = unique_instance_gui_msg__;
                % this.getUniqueWinHandle();
                if nargin == 0 || ~flag_skip_init
                    this.init();
                    this.openGUI();
                else
                    this.openGUI(flag_skip_init);
                end                
            end
            
            drawnow
        end
    end
    
    %% METHODS INIT
    % ==================================================================================================================================================
    methods
        function init(this)
        end
        
        function openGUI(this, flag_skip_init)
            % Main Window ----------------------------------------------------------------------------------------------
            
            % If there is still an old logging wondow still open, close it
            old_win = this.getUniqueWinHandle();
            if nargin > 1 && flag_skip_init && old_win.isvalid
                this.win = old_win;
            else
                if ~isempty(old_win)
                    delete(old_win);
                end
                
                win = figure( 'Name', 'goGPS log', ...
                    'Visible', 'off', ...
                    'DockControls', 'off', ...
                    'MenuBar', 'none', ...
                    'ToolBar', 'none', ...
                    'NumberTitle', 'off', ...
                    'Resize', 'on');
                win.Position([3,4]) = [500 640]; % Set dimensions
                win.UserData.name = this.WIN_NAME;
                
                this.win = win;
                
                if isunix && not(ismac())
                    % top right
                    % win.Position(1) = round((win.Parent.ScreenSize(3) - win.Position(3)));
                    % win.Position(2) = round((win.Parent.ScreenSize(4) - win.Position(4)));
                    % centered
                    win.Position(1) = 0; %round((win.Parent.ScreenSize(3) - win.Position(3)) / 2);
                    win.Position(2) = round((win.Parent.ScreenSize(4) - win.Position(4)) / 2);
                else
                    % top right
                    % win.OuterPosition(1) = round((win.Parent.ScreenSize(3) - win.OuterPosition(3)));
                    % win.OuterPosition(2) = round((win.Parent.ScreenSize(4) - win.OuterPosition(4)));
                    % centered
                    win.OuterPosition(1) = 5; %round((win.Parent.ScreenSize(3) - win.OuterPosition(3)) / 2);
                    win.OuterPosition(2) = round((win.Parent.ScreenSize(4) - win.OuterPosition(4)) / 2);
                end
                
                try
                    main_vb = uix.VBox('Parent', win, ...
                        'Padding', 5, ...
                        'BackgroundColor', Core_UI.DARKER_GREY_BG);
                catch
                    log = Core.getLogger;
                    log.setOutMode(1,[],0); % to plot a Warning I need to disable GUI and enable
                    log.addError('Please install GUI Layout Toolbox (https://it.mathworks.com/matlabcentral/fileexchange/47982-gui-layout-toolbox)');
                    open('GUI Layout Toolbox 2.3.4.mltbx');
                    log.newLine();
                    log.addWarning('After installation re-run goGPS');
                    close(win);
                    return;
                end
                top_bh = uix.HBox('Parent', main_vb);
                
                logo_GUI_Msg.BG_COLOR = Core_UI.DARK_GREY_BG;
                left_tbv = uix.VBox('Parent', top_bh, ...
                    'BackgroundColor', logo_GUI_Msg.BG_COLOR, ...
                    'Padding', 5);
                
                % Logo/title box -------------------------------------------------------------------------------------------
                
                logo_g = uix.Grid('Parent', left_tbv, ...
                    'Padding', 5, ...
                    'BackgroundColor', logo_GUI_Msg.BG_COLOR);
                
                logo_ax = axes( 'Parent', logo_g);
                logo_g.Widths = 64;
                logo_g.Heights = 64;
                [logo, transparency] = Core_UI.getLogo();
                logo(repmat(sum(logo,3) == 0,1,1,3)) = 0;
                logo = logo - 20;
                image(logo_ax, ones(size(logo)), 'AlphaData', transparency);
                logo_ax.XTickLabel = [];
                logo_ax.YTickLabel = [];
                axis off;
                
                Core_UI.insertEmpty(left_tbv, logo_GUI_Msg.BG_COLOR);
                left_tbv.Heights = [82 -1];
                
                % Title Panel -----------------------------------------------------------------------------------------------
                right_tvb = uix.VBox('Parent', top_bh, ...
                    'Padding', 5, ...
                    'BackgroundColor', logo_GUI_Msg.BG_COLOR);
                
                top_bh.Widths = [106 -1];
                
                title = uix.HBox('Parent', right_tvb, ...
                    'BackgroundColor', logo_GUI_Msg.BG_COLOR);
                
                txt = this.insertBoldText(title, 'goGPS', 12, Core_UI.LBLUE, 'left');
                txt.BackgroundColor = logo_GUI_Msg.BG_COLOR;
                title_l = uix.VBox('Parent', title, 'BackgroundColor', GUI_Msg.BG_COLOR);
                title.Widths = [60 -1];
                Core_UI.insertEmpty(title_l, logo_GUI_Msg.BG_COLOR)
                txt = this.insertBoldText(title_l, ['- software V' Core.GO_GPS_VERSION], 9, [], 'left');
                txt.BackgroundColor = logo_GUI_Msg.BG_COLOR;
                title_l.Heights = [2, -1];
                
                % Disclaimer Panel -----------------------------------------------------------------------------------------------
                Core_UI.insertEmpty(right_tvb, logo_GUI_Msg.BG_COLOR)
                txt = Core_UI.insertText(right_tvb, {'A GNSS processing software powered by GReD'}, 9, GUI_Msg.BG_COLOR, [], 'left');
                txt.BackgroundColor = logo_GUI_Msg.BG_COLOR;
                right_tvb.Heights = [25 3 -1];
                
                % Logging Panel --------------------------------------------------------------------------------------------------
                log_container = uix.VBox('Parent', main_vb, 'Padding', 5, 'BackgroundColor', GUI_Msg.BG_COLOR);
                [j_edit_box, h_log_panel] = Core_UI.insertLog(log_container);
                this.win.UserData.jedt = j_edit_box;
                this.jedt = j_edit_box;
                
                this.clear();
                
                % Manage dimension -------------------------------------------------------------------------------------------
                
                main_vb.Heights = [84 -1];
            end
            this.win.Visible = 'on';
        end
        
        function close(this)
            if ~isempty(this.win) && ishandle(this.win)
                close(this.win);
            end
        end
    end
    %% METHODS INSERT
    % ==================================================================================================================================================
    methods (Static)
        function txt = insertBoldText(parent, title, font_size, color, alignment)
            if nargin < 4 || isempty(color)
                color = Core_UI.WHITE;
            end
            if nargin < 5 || isempty(alignment)
                alignment = 'center';
            end
            txt = uicontrol('Parent', parent, ...
                'Style', 'Text', ...
                'String', title, ...
                'ForegroundColor', color, ...
                'HorizontalAlignment', alignment, ...
                'FontSize', Core_UI.getFontSize(font_size), ...
                'FontWeight', 'bold', ...
                'BackgroundColor', GUI_Msg.BG_COLOR);
        end
    end
    %% METHODS SETTERS
    % ==================================================================================================================================================
    methods
        function addMessage(this, text, type)
            % Add a message to the logger
            % 
            % INPUT
            %   text    text in HTML format
            %   type    'm'     marked message
            %           'w'     warning message
            %           'e'     error message
            %           otherwise normal
            %
            % SYNTAX
            %   this.addHTML(text, type)
            
            this.bringOnTop();
            if nargin < 3 || isempty(type)
                type = 'n';
            end
            
            Core_UI.guiAddMessage(this.jedt, text, type);
        end
        
        function addHTML(this, html_txt)
            % Add a message to the logger
            % 
            % INPUT
            %   html_txt    text in HTML format
            %
            % SYNTAX
            %   this.addHTML(html_txt)

            this.bringOnTop();
            Core_UI.guiAddHTML(this.jedt, html_txt);
        end
        
        function clear(this)
            % Add a message to the logger
            % 
            % SYNTAX
            %   this.clear()
            
            this.bringOnTop();
            Core_UI.guiClearLog(this.jedt);
            Core_UI.guiAddMessage(this.jedt, ['<p>' GPS_Time.now.toString('yyyy/mm/dd HH:MM:SS') '</p>']);
            Core_UI.guiAddMessage(this.jedt, ['<p><b>Welcome to goGPS!</b></p>for any problem contact us at <a color="' rgb2hex(Core_UI.LBLUE) '" source="http://bit.ly/goGPS">http://bit.ly/goGPS</a>'], 'm');            
            this.addHTML(['<font color=gray face="Courier">' repmat('&mdash;', 1, 56) ' </font>']);
        end
    end
    
    %% METHODS STATIC 
    % ==================================================================================================================================================
    methods (Static)
        function pos = getPosition()
            pos = GUI_Msg.getInstance(true).win.Position;
        end
    end
    %% METHODS EVENTS
    % ==================================================================================================================================================
    methods (Access = public)         
    end
end
