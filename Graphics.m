
% Created by Jorge M. Cruz-Duarte, feb-2018
% Contact: jorge.cruz@ugto.mx
%
% Example:
% >> q = Graphics('Hello_World'); % Create the object
%
% >> plot(rand(1,3), rand(1,3),'b'), hold on,
% >> plot(rand(1,3), rand(1,3),'r'),
% >> xlabel('x-axis'); ylabel('y-axis');
% >> legend('Data 1','Data 2');
%
% >> setup(q); % Quick setup
% >> save(q); % Print graphic
%
classdef Graphics
    properties
        objID
        fileName
    end
    methods
        %% Construct an object and initialise it
        function obj = Graphics(namestring)
            if nargin > 0
                if ischar(namestring)
                    % Crate the object
                    obj.fileName = namestring;
                    obj.objID = figure('color','w','name',namestring);
                    % Set default properties
                    setup(obj);
                else
                    error('Name must be char!')
                end
            else
                [~,namestring] = fileparts(tempname);
                obj = Graphics(namestring);
            end
        end
        
        %% Destroy an object
        function kill(obj)
            if nargin > 0
                close(obj.objID);
            else
                error('An object is needed!!!');
            end
        end        
        
        %% Save the graphic as an eps figure
        function save(obj,otherName,extFile)
            delim  = '_';
            if nargin < 3
                extFile = '.eps';
                if nargin < 2
                    otherName = '';
                    delim = '';
                end
            end
            if strcmp(otherName,'')
               name = obj.fileName;
            else
                name = [obj.fileName,delim,otherName];
            end
            
            switch extFile
                case '.jpg'
                    print(obj.objID,[name,extFile],'-r333','-djpeg','-noui');
                case '.png'
                    print(obj.objID,[name,extFile],'-r333','-dpng','-noui');
                otherwise
                    print(obj.objID,[name,extFile],'-painters','-depsc','-tiff','-noui');
            end
        end
        
        %% Setting some properties of the graphics
        function setsize(obj,columns,aspectRatio)
            if nargin > 0
                % Get paper size
                paperSize = [21.59,27.94]; % in cm %get(obj.objID,'PaperSize');
                % Set margins
                margin = 2; % in cm
                columnsep = 0.35; % in cm
                % textFullWidth
                textFullWidth = paperSize(2) - 2*margin;
                
                warningNargin1 = 'Nothing was set! Introduce an object';
                
                if nargin > 1
                    % Read Columns
                    if columns > 1 %&& mod(columns,1) == 0
                        lineWidth = textFullWidth/columns - columnsep/2;
                    else
                        lineWidth = textFullWidth;
                    end
                    
                    % Read aspectRatio
                    errorReadAR.message = ['Aspect ratio must be char, ',...
                        'number or 2D vector. For example: ',...
                        '''4:2 (default)'', ','0.25, or [12,4]'];
                    errorReadAR.identifier = 'Graphics:aspectRatioRead';
                    
                    if ischar(aspectRatio)
                        asprs = split(aspectRatio,':');
                        if numel(asprs) ~= 2, error(errorReadAR); end
                        
                        alpha = 1/eval(sprintf('%s/%s',asprs));
                    elseif isnumeric(aspectRatio)
                        switch numel(aspectRatio)
                            case 1
                                alpha = aspectRatio;
                            case 2
                                alpha = aspectRatio(2)/aspectRatio(1);
                            otherwise
                                error(errorReadAR);
                        end
                    else
                        error(errorReadAR);
                    end
                    
                    % Set figure size
                    graphicSize = lineWidth*[1,alpha]; % width, height in cm
                    %                     set(obj.objID,'Menubar','none');
                    set(obj.objID,'PaperUnits','centimeters');
                    set(obj.objID,'PaperSize',graphicSize);
                    set(obj.objID,'PaperPosition',[0 0 graphicSize]);
                    set(obj.objID,'Units','centimeters',...
                        'Position',[0 0 graphicSize]);
                    movegui(obj.objID, 'center');
                    
                else % Default set size: columns = 1, aspectRatio = 2:1
                    setsize(obj,1,[2,1]);
                end
            else
                warning(warningNargin1)
            end
        end
        
        %% Setting font size
        function setfont(obj,fontSize)
            if nargin > 0
                if nargin > 1
                    
                    % Check fontSize
                    errorReadFS.message = ['Font size must be number',...
                        ' in pt. For example: 10.5 or 12 (default)'];
                    errorReadFS.identifier = 'Graphics:fontSizeRead';
                    if ~(isnumeric(fontSize) && numel(fontSize) == 1)
                        error(errorReadFS);
                    end
                    
                    % Set fonts
                    if ~isempty(obj.objID.CurrentAxes)
                        
                        % Find all axes
                        allAxes = findall(obj.objID,'type','axes');
                        
                        for ia = 1 : numel(allAxes)
                            ax = allAxes(ia);
                            
                            set(ax,'FontSize',fontSize);
                            set(ax,'TickLabelInterpreter','latex');
                            set(ax.XLabel,'Interpreter','latex');
                            set(ax.YLabel,'Interpreter','latex');
                            set(ax.ZLabel,'Interpreter','latex');
                            
                            if ~isempty(ax.Legend)
                                set(ax.Legend,'FontSize',fontSize);
                                set(ax.Legend,'Interpreter','latex');
                            end
                        end
                    end
                else
                    setfont(obj,12)
                end
            else
                warning(warningNargin1);
            end
        end
        
        %% Set line
        function setline(obj,lineWidth)
            
            if nargin > 1
                errorReadLW.message = ['Line width must be number',...
                    ' in pt. For example: 0.5 or 1 (default)'];
                errorReadLW.identifier = 'Graphics:fontSizeRead';
                if ~(isnumeric(lineWidth) && numel(lineWidth) == 1)
                    error(errorReadLW);
                end
                
                if ~isempty(obj.objID.Children)
                    
                    % Find all axes
                    allAxes = findall(obj.objID,'type','axes');
                    
                    for ia = 1 : numel(allAxes)
                        ax = allAxes(ia);
                        set(ax,'LineWidth',lineWidth);
                    end
                end
            else
                setline(obj,1)
            end
        end
        
        %% quick setup
        function setup(obj)
            setsize(obj);
            setfont(obj);
            setline(obj);
            
            TextEPS = [
                1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
                1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
                1 1 1 1 1 0 0 1 1 1 1 1 0 0 0 0
                1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0
                1 0 0 0 0 0 0 1 0 0 0 1 0 0 0 0
                1 1 1 1 1 0 0 1 1 1 1 1 0 0 0 0
                0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0
                0 0 0 1 0 0 0 1 0 0 1 1 1 1 1 1
                0 0 0 1 0 0 0 1 0 0 1 0 0 0 0 0
                0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0
                1 0 0 1 0 0 1 0 0 0 1 1 1 1 1 1
                0 1 0 1 0 1 0 0 0 0 0 0 0 0 0 1
                0 0 1 1 1 0 0 0 0 0 0 0 0 0 0 1
                1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1
                ];
            TextEPS(TextEPS == 0) = nan;
            
            TextPNG = [
                1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0
                1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
                1 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0
                1 1 1 1 1 0 1 0 0 0 0 0 0 0 0 0
                1 0 0 0 0 0 1 1 0 0 0 1 0 0 0 0
                1 0 0 0 0 0 1 0 1 0 0 1 0 0 0 0
                1 0 0 0 0 0 1 0 0 1 0 1 0 0 0 0
                1 0 0 0 0 0 1 0 0 0 1 1 0 0 0 0
                0 0 0 0 0 0 1 0 0 0 0 1 0 0 0 0
                0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 1 0 0 0 0 0 0 1 1 1 1 1 1
                0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0
                1 0 0 1 0 0 1 0 0 0 1 0 0 0 0 0
                0 1 0 1 0 1 0 0 0 0 1 0 0 1 1 1
                0 0 1 1 1 0 0 0 0 0 1 0 0 0 0 1
                1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1
                ];
            TextPNG(TextPNG == 0) = nan;
            
            TextJPG = [
                0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0
                0 0 0 1 0 0 1 1 1 1 1 0 0 0 0 0
                0 0 0 1 0 0 1 0 0 0 1 0 0 0 0 0
                0 0 0 1 0 0 1 0 0 0 1 0 0 0 0 0
                1 1 1 1 0 0 1 1 1 1 1 0 0 0 0 0
                0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
                0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0
                0 0 0 1 0 0 1 0 0 0 0 0 0 0 0 0
                0 0 0 1 0 0 0 0 0 0 1 1 1 1 1 1
                0 0 0 1 0 0 0 0 0 0 1 0 0 0 0 0
                1 0 0 1 0 0 1 0 0 0 1 0 0 0 0 0
                0 1 0 1 0 1 0 0 0 0 1 0 0 1 1 1
                0 0 1 1 1 0 0 0 0 0 1 0 0 0 0 1
                1 1 1 1 1 1 1 0 0 0 1 1 1 1 1 1
                ];
            TextJPG(TextJPG == 0) = nan;
            
            CDataEPS = nan(16,16,3);
            CDataEPS(:,:,1) = TextEPS;
            
            CDataPNG = nan(16,16,3);
            CDataPNG(:,:,1) = TextPNG;
            CDataPNG(:,:,3) = TextPNG;
            
            CDataJPG = nan(16,16,3);
            CDataJPG(:,:,1) = 0.2*TextJPG;
            CDataJPG(:,:,3) = TextJPG;
            
            % Put a uibutton for print
            if isempty(findall(obj.objID,'Tag','.eps'))
                uipushtool(findall(obj.objID,'Type','uitoolbar'),...
                    'CData',CDataEPS,...
                    'TooltipString','Save to EPS', 'Tag', '.eps', ...
                    'Separator','on','HandleVisibility','off',...
                    'ClickedCallback',@pushbutton1_Callback);
            end
            
            if isempty(findall(obj.objID,'Tag','.png'))
                uipushtool(findall(obj.objID,'Type','uitoolbar'),...
                    'CData',CDataPNG,...
                    'TooltipString','Save to PNG', 'Tag', '.png', ...
                    'Separator','on','HandleVisibility','off',...
                    'ClickedCallback',@pushbutton1_Callback);
            end
            
            if isempty(findall(obj.objID,'Tag','.jpg'))
                uipushtool(findall(obj.objID,'Type','uitoolbar'),...
                    'CData',CDataJPG,...
                    'TooltipString','Save to JPG', 'Tag', '.jpg', ...
                    'Separator','on','HandleVisibility','off',...
                    'ClickedCallback',@pushbutton1_Callback);
            end
            
            %             ButtonH = uicontrol('Parent',obj.objID,'Style','Pushbutton',...
            %                 'String','Print','Units','Normalized',...
            %                 'Position',[0.88 0.9 0.11 0.07],'Visible','on',...
            %                 'Callback',@pushbutton1_Callback);
            
            %             function Pushing()
            %                 [~,RandomEnding] = fileparts(tempname);
            %                 ButtonH.Visible = 'off';
            %                 save(gcf,RandomEnding(end-11:end));
            %                 ButtonH.Visible = 'on';
            %             end
            function pushbutton1_Callback(hObject, eventdata, handles)
                fprintf('Saving...\n');
                % hObject.Tag = .eps | .jpg | .png
                save(obj,'',hObject.Tag);
                fprintf('Saved as ''%s%s''\n',obj.fileName,hObject.Tag);
            end
        end
        
        % set all up
        function setall(obj,columns,aspectratio,font,line)
            setup(obj);
            setsize(obj,columns,aspectratio);
            setfont(obj,font);
            setline(obj,line);
        end
    end
end % class
