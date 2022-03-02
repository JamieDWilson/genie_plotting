%% plot_genie_lonlat
% -------------------------------------------------------------------------
% Description: Quick interactive plots of lon/lat fields from genie output.
%
% Inputs: 
%   - output_dir: name or path of genie output file, e.g., 'fields_biogem_3d.nc' or 'experiment_dir/biogem/fields_biogem_3d.nc'
%   - var_name: name of variable in netCDF to plot, e.g., 'ocn_PO4'
%   - year: timeslice to plot, e.g., 9999.5
%   - depth: depth level to plot (1 for surface, >1 for interior)
%   - data_scale: scale output by orders of magnitude, e.g., 1e-6 (micro)
%   - optional #1: name or path to 2nd output file to plot as difference
%   - optional #2: timeslice of 2nd output to plot as difference
%
% Outputs:
%   - object handle, used to interactively change the plot and/or call subfunctions
%
% User options:
%   - autoplot: redraws figure when initiating object and changing parameters (logical)
%   - cmin: minimum value for colorscale (float)
%   - cmax: maximum value for colorscale (float)
%   - c_nlevels: number of colour intervals (integer)
%   - colormap: name of colormap (string)
%   - colorbar: plot a colourbar (logical)
%   - colorbar_text: colorbar text (string)
%   - title_text: text for plot title (string)
%   - lon_label: have a longitude label (logical)
%   - lat_label: have a latitude label (logical)
%   - lat_tick_labels: latitude tick labels to display (float)
%   - lon_tick_labels: longitude tick labels to display (float)
%   - lon_origin: longitude of origin for plot 
%   - overlay_data: point data to overlay as a mx3 array with cols of lon, lat, data (float)
%   - overlay_point_size: size of overlay points (integer)
%
% Subfunctions:
%   - calc_overlay_grid_val: calculates lon/lat/value of overlay_data locations on genie grid
%   - plot_zonal_mean: plots a zonal mean of data. Line properties can be passed as name-value pairs to customise plot
%   - save_bitmap: saves plot to a .png file
%   - save_vector: saves plot to a .svg file
%
% Usage / Examples:
%
%   --------------------- INTERACTIVE PLOTTING ----------------------------
%
%   - with autoplot=true the plot is automatically drawn upon calling
%   plot_fields_lonlat() and re-drawn when changing parameters so that you
%   can interactively alter the plot. If more than plot is open, the object
%   the plot will be re-drawn on its respective original figure window.
%
%   -> plot basic PO4 lon/lat field and adjust colourscale interactively
%   fig=plot_genie_lonlat_alt ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 )
%   fig.cmin=0;
%   fig.cmax=2;
%   fig.c_n_levels=10;
%
%   -> plot difference with another experiment
%   fig=plot_genie_lonlat_alt ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 , 'fields_biogem_3d_alt.nc' , 9999.5 )
%
%   -> plot field, overlay datapoints, and retrieve values from genie grid
%   fig=plot_genie_lonlat_alt ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 )
%   fig.overlay_data=overlay_data_array;
%   grid_data=fig.calc_overlay_grid_val;
%
%   -> plot field and zonal mean
%   fig=plot_genie_lonlat_alt ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 )
%   figure;
%   fig.plot_zonal_mean;
%
%   --------------------- STATIC PLOTTNG & SUBPLOTs -----------------------
%
%   - setting autoplot=false allows different usage. Parameter values can
%   be set and the plot drawn finally using the .plot function. This allows
%   plots to be drawn as subplots. This draws to the currently
%   selected figure window.
%
%   -> plot single final figure, e.g., in a script
%   fig=plot_genie_lonlat_alt ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 )
%   fig.cmin=0;
%   fig.cmax=2;
%   fig.c_n_levels=10;
%   fig.plot;
%   
%   -> plot as a subplot
%   figure
%   subplot(2,1,1)
%   fig=plot_genie_lonlat_alt ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 )
%   fig.cmin=0;
%   fig.cmax=2;
%   fig.c_n_levels=10;
%   fig.plot;
%   subplot(2,1,2)
%   fig=plot_genie_lonlat_alt ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 )
%   fig.cmin=0;
%   fig.cmax=2;
%   fig.c_n_levels=10;
%   fig.plot;
%   
% Notes:
%   - calling the main function will produce the basic plot
%   - everytime you change user-definable options, the plot will be redrawn
%   - adding an additional optional output will display output #1 - output #2
%   - if a variable name, year or depth is not available, a list of available ones will be printed to screen
%   - the colour scale is by default fitted to the data
%
% Known Issues:
%   - overlay data might plot off the map if lon/lat are close to the edges
% -------------------------------------------------------------------------



classdef plot_genie_lonlat < handle
    
    methods
        
        function obj=plot_genie_lonlat ( output_file , var_name , year , depth , data_scale , varargin )
            
            % set class inputs to variables
            obj.output_file = output_file; 
            obj.var_name = var_name;   
            obj.year = year;   
            obj.depth = depth;
            obj.data_scale = data_scale;
            obj.opt_args = varargin;
            
            % load in plot data
            initialise_plot_data(obj);
            
            % initial plot
            if obj.autoplot
                plot(obj);
            end
            
        end
        
    end
    
    properties(SetObservable) 
        
        %---------------------------------%
        % --- USER DEFINABLE DEFAULTS --- %
        %---------------------------------%
        
        % autoplot when calling function and changing parameters?
        % if false, use .plot to draw figure
        % (plots to the current figure/axes handle)
        autoplot=true;
        
        %---------------------------------%
        % ------ DEFAULT PARAMETERS ----- %
        %---------------------------------%
        
         %--------- COLORSCALE ----------%
        % minimum colour ([] = min of data)
        cmin=[];
        % maximum colour ([] = max of data)
        cmax=[];
        % number of colour degrees
        c_nlevels=200;
        % colormap ('auto' = parula for single dataset, RdBu for difference)
        colormap='auto';
        
        %---------- COLORBAR ----------%
        % colorbar?
        colorbar=true;
        % colorbar text ('auto' = uses variable name and units from netcdf)
        colorbar_text = {'auto'};
        
        %---------- DISPLAY -----------%
        % plot title ('auto' = automatic using file name, variable name, depth and year)
        title_text = {'auto'};
        % longitude label?
        lon_label=true;
        % latitude label?
        lat_label=true;
        % latitude tick labels (deg N)
        lat_tick_label=[-60 -30 -0 30 60];
        % longitude tick labels (deg E)
        lon_tick_label=[-360:60:360];
        % longitude origin (in deg E, [] = origin in file)
        lon_origin=[];
        
        %-------- OVERLAY DATA --------%
        % point dataset to overlay (default - no overlay data)
        overlay_data=[];
        % size of overlay datapoints
        overlay_point_size=60;
         
    end
    
    properties(GetAccess=public,SetAccess=public)

        % ---- placeholders ----
        output_file = 1; 
        var_name  = 1;   
        year = 1;  
        depth  = 1;
        data_scale  = 1;
        opt_args  = 1;
        data_n = 1;
        data_2_flag = 1;
        output_dirs = {};
        lat = 1;
        lon = 1;
        z = 1;
        time = {};
        units = 1;
        longname = 1;
        data = 1;
        data2 = 1;
        c = 1;
        coastlines = 1;
        n_lon = 1;
        n_lat = 1;
        lon_to_i = 1;
        lat_to_j = 1;
        si_om = 1;
        lon_ticks = 1;
        lat_ticks = 1;
        lat_to_y = 1;
        lon_to_x = 1;
        overlay_data_x = 1;
        overlay_data_y = 1;
        fig=1;
        ax=1;
        im=1;
        lon_origin_data=1;
        
    end
        
        

    methods(Access = public)
        
        function initialise_plot_data(obj)
            
            % check 2nd data set details are entered correctly
            if ~isempty(obj.opt_args) & size(obj.opt_args,2)~=2
                disp('Wrong number of arguments for 2nd timeslice')
                disp('Last 2 inputs need to be:') 
                disp('output_directory , year')
                return
            end

            % assign 2nd data set info
            obj.data_n=1;
            obj.data_2_flag=false;
            obj.output_dirs{1}=obj.output_file;
            if ~isempty(obj.opt_args)
                obj.output_dirs{2}=obj.opt_args{1};
                obj.year=[obj.year; obj.opt_args{2}];
                obj.data_n=2;
                obj.data_2_flag=true;
            end

            % grid variables
            obj.lat=ncread(obj.output_dirs{1},'lat');
            obj.lon=ncread(obj.output_dirs{1},'lon');
            obj.z=ncread(obj.output_dirs{1},'zt');
            obj.n_lon=numel(obj.lon);
            obj.n_lat=numel(obj.lat);
            obj.lon_origin_data=min(obj.lon)-unique(diff(obj.lon))/2;
            %area=ncread(obj.output_dirs{1},'grid_A');
            %dD=ncread(obj.output_dirs{1},'grid_dD');
            %vol=area.*dD;
            %obj.kg=1027.649.*vol;

            % time variables
            for n=1:obj.data_n
                obj.time{n,1}=ncread(obj.output_dirs{n},'time');
            end

            % check year exists
            for n=1:obj.data_n
                if ~ismember(obj.year(n),obj.time{n})
                    disp('Year selected is not available. Choose from:')
                    disp(num2str(obj.time{n}))
                    return
                end
            end

            % find output index corresponding to year
            for n=1:obj.data_n
                obj.year(n)=find(obj.year(n)==obj.time{n}); % convert to index
            end

            % check depth selection
            if obj.depth>numel(obj.z)
                disp('Depth level selected is larger than available')
                disp(['Maximum levels available: ' num2str(numel(obj.z))])
                return
            end

            % check variable exists & get metadata
            for nn=1:obj.data_n

                netcdf_info=ncinfo(obj.output_dirs{nn});
                netcdf_info=netcdf_info.Variables;

                var_nc_ind=find(strcmp({netcdf_info.Name},obj.var_name));
                if isempty(var_nc_ind)
                    disp('Variable name not available')
                    disp(['Choose from: '])
                    for n=1:size(netcdf_info,2)
                        disp(netcdf_info(n).Name)
                    end
                    return
                end

                if any(strcmp({netcdf_info(var_nc_ind).Attributes.Name},'units'))
                    ind=find(strcmp({netcdf_info(var_nc_ind).Attributes.Name},'units'));
                    obj.units=netcdf_info(var_nc_ind).Attributes(ind).Value;
                else
                    obj.units='?';
                    warning('No units found for variable')
                end
                if any(strcmp({netcdf_info(var_nc_ind).Attributes.Name},'long_name'))
                    ind=find(strcmp({netcdf_info(var_nc_ind).Attributes.Name},'long_name'));
                    obj.longname=netcdf_info(var_nc_ind).Attributes(ind).Value;
                else
                    obj.longname='?';
                    warning('No longname found for variable')
                end

            end

            % get data
            obj.data=ncread(obj.output_dirs{1},obj.var_name);
            if obj.data_2_flag
                obj.data2=ncread(obj.output_dirs{2},obj.var_name);
            end

            % set missing values to NaNs
            obj.data(obj.data>1e30)=NaN;

            % slice data
            netcdf_info=ncinfo(obj.output_dirs{1});
            var_nc_ind=find(strcmp({netcdf_info.Variables.Name},obj.var_name));
            % check whether 3D 
            if any(strcmp({netcdf_info.Variables(var_nc_ind).Dimensions.Name},'zt'))
                obj.data=squeeze(obj.data(:,:,obj.depth,obj.year(1)));
                if obj.data_2_flag
                    obj.data2=squeeze(obj.data2(:,:,obj.depth,obj.year(2)));
                    obj.data=obj.data-obj.data2;
                end
            else % assume 2D
                obj.data=squeeze(obj.data(:,:,obj.year(1)));
                if obj.data_2_flag
                    obj.data2=squeeze(obj.data2(:,:,obj.year(2)));
                    obj.data=obj.data-obj.data2;
                end
            end

            % scale data
            obj.data=obj.data./obj.data_scale;
            
            % orient data lon (x) lat (y)
            obj.data=rot90(obj.data);
            
            % colour scale
            if obj.data_2_flag
                maxmin=max([abs(min(min(obj.data))) abs(max(max(obj.data)))]);
                obj.cmin=-maxmin;
                obj.cmax=maxmin;
            else
                obj.cmin=min(min(obj.data));
                obj.cmax=max(max(obj.data));
            end
            
            % orders of magnitude symbols
            obj.si_orders_of_mag;
            
            % auto title text
            if strmatch(obj.title_text','auto')
                if obj.data_2_flag
                    obj.title_text={ [obj.longname ' at ' num2str(round(obj.z(obj.depth),2)) ' m'] , ['(' obj.output_dirs{1} ' @ ' num2str(obj.time{1}(obj.year(1))) ' years) - (' obj.output_dirs{2} ' @ ' num2str(obj.time{2}(obj.year(2))) ' years)'] };
                else
                    obj.title_text={ [obj.longname ' at ' num2str(round(obj.z(obj.depth),2)) ' m'] , ['(' obj.output_dirs{1} ' @ ' num2str(obj.time{1}(obj.year(1))) ' years)'] };
                end
            end
            
            % auto colorbar text
            if strmatch(obj.colorbar_text,'auto')
                obj.colorbar_text=[obj.longname ' (' obj.si_om obj.units ')'];
            end
            
            % figure handle
            if obj.autoplot
                obj.fig=figure;
            end
            
            % everything is now initialised
            % add listeners to monitor changes in certain parameters
            % from now on changing these parameters will redraw the plot
            addlistener(obj,'c_nlevels','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'lat_tick_label','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'lon_tick_label','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'cmin','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'cmax','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'lon_origin','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'overlay_data','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'title_text','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'colorbar_text','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'colorbar','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'lon_label','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'lat_label','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'overlay_point_size','PostSet',@obj.handlePropertyEvents);
            addlistener(obj,'colormap','PostSet',@obj.handlePropertyEvents);
            
        end
        
        function [] = plot(obj)
               
            % change origin
            [plot_data,plot_lon]=obj.change_origin;
            
            % make coastlines
            obj.make_coastlines(plot_data);
            
            % make colormap
            obj.make_colormap;
            
            % point at relevant figure handle
            if obj.autoplot
                figure(obj.fig);
            end

            % plot data, NaNs not plotted 
            clf
            imAlpha=ones(size(plot_data));
            imAlpha(isnan(plot_data))=0;
            obj.im=imagesc([1:36],[1:36],plot_data,'AlphaData',imAlpha);
            
            % get axes
            if ~obj.autoplot
                obj.ax=gca;
            else
                obj.ax=obj.fig.CurrentAxes;
            end
            
            % lon/lat axes formatting
            obj.lonlat_to_xy(plot_lon,obj.lat);
            obj.plot_ticks;
            set(obj.ax,'YTick',obj.lat_ticks,'YTickLabel',num2cell(-obj.lat_tick_label),'XTick',obj.lon_ticks,'XTickLabel',num2cell(obj.lon_tick_label)) % set custom ticks/labels
            if obj.lon_label
                xlabel('Longitude')
            end
            if obj.lat_label
                ylabel('Latitude')
            end
            
            % colorbar
            if obj.colorbar
                C=colorbar;
                ylabel(C,obj.colorbar_text);
            end
            
            % plot formating
            axis normal
            pbaspect([2 1 1]) % rectangular aspect
            hold on
    
            % plot coastlines
            for i=1:numel(obj.coastlines(:,1))
                plot(obj.coastlines(i,1:2)-1.5,obj.coastlines(i,3:4)-1.5,'k')
                hold on
            end
            
            % set background (and land) color
            set(obj.ax, 'color', [0.8 0.8 0.8]);
            
            % colormaps
            colormap(obj.c);
            caxis(obj.ax,[obj.cmin obj.cmax]);
            
            % title
            title(obj.title_text,'Interpreter','None');
            
            % overlay data
            if ~isempty(obj.overlay_data)
                
                obj.process_overlay_data(plot_lon);

                if isnan(obj.overlay_data(:,3))
                    scatter(obj.overlay_data_x,obj.overlay_data_y,obj.overlay_point_size,ones(numel(obj.overlay_data(:,3)),3)*0.6,'filled','MarkerEdgeColor','k');
                else
                    scatter(obj.overlay_data_x,obj.overlay_data_y,obj.overlay_point_size,obj.overlay_data(:,3),'filled','MarkerEdgeColor','k');
                end
            end
            
            
        end
    
        function [] = make_colormap(obj)
            
            if strmatch(obj.colormap,'auto')
                if obj.data_2_flag
                  obj.c=obj.make_cmap('RdBu',obj.c_nlevels);
                else
                  obj.c=obj.make_cmap('parula',obj.c_nlevels);
                end
            else
                obj.c=obj.make_cmap(obj.colormap,obj.c_nlevels);
            end
            
        end
        
        function [] = make_coastlines(obj,data)
           
            
            grid_mask=zeros(obj.n_lon,obj.n_lat);
            grid_mask(~isnan(data))=1;
            grid_mask = [grid_mask(:,1) grid_mask grid_mask(:,end)];
            grid_mask = [grid_mask(1,:) ;grid_mask; grid_mask(end,:)];

            dx=diff(grid_mask,1,2); 
            dy=diff(grid_mask,1,1); 
            [jv iv]=find(dx~=0);
            [jh ih]=find(dy~=0);
            vertedge=[iv+1 iv+1 jv jv+1];
            horzedge=[ih ih+1 jh+1 jh+1];
            obj.coastlines=[vertedge;horzedge];
            
        end
        
        function [] = lonlat_to_xy(obj,lon,lat)
            
            obj.lon_to_x=griddedInterpolant(lon,[1:1:obj.n_lon]);
            obj.lat_to_y=griddedInterpolant(lat,[obj.n_lat:-1:1]);
            
        end
        
        function [] = si_orders_of_mag(obj)
            
            % scale colorbar units
            switch obj.data_scale
                case 1e3
                    obj.si_om='k';
                case 1e6
                    obj.si_om='M';
                case 1e9
                    obj.si_om='G';
                case 1e12
                    obj.si_om='T';
                case 1e15
                    obj.si_om='P';
                case 1e-3
                    obj.si_om='m';
                case 1e-6
                    obj.si_om='\mu';
                case 1e-9
                    obj.si_om='n';
                case 1e-12
                    obj.si_om='p';
                otherwise
                    obj.si_om='';
            end
        end
        
        function [] = plot_ticks(obj)
            
            obj.lat_ticks=fliplr(obj.lat_to_y(obj.lat_tick_label));
            obj.lon_ticks=obj.lon_to_x(obj.lon_tick_label);

        end
        
        function [] = process_overlay_data(obj,lon)
            
            % change longitude values 
            obj.overlay_data_x=obj.overlay_data(:,1);
            obj.overlay_data_x(obj.overlay_data_x<min(lon))=obj.overlay_data_x(obj.overlay_data_x<min(lon))+360;
            obj.overlay_data_x(obj.overlay_data_x>max(lon))=obj.overlay_data_x(obj.overlay_data_x>max(lon))-360;
            
            obj.overlay_data_x=obj.lon_to_x(obj.overlay_data_x);
            obj.overlay_data_y=obj.lat_to_y(obj.overlay_data(:,2));
            obj.overlay_data(:,3)=obj.overlay_data(:,3)./obj.data_scale;
            
        end

       function [overlay_data_grid] = calc_overlay_grid_val(obj)
           
           overlay_data_grid=zeros(size(obj.overlay_data));
           
           for n=1:numel(obj.overlay_data(:,1))
               
                [~,lon_id]=min(abs(obj.lon-obj.overlay_data(n,1)));
                [~,lat_id]=min(abs(obj.lat-obj.overlay_data(n,2)));
                
                overlay_data_grid(n,1)=obj.lon(lon_id);
                overlay_data_grid(n,2)=obj.lat(lat_id);
                
                x=obj.lon_to_x(obj.lon(lon_id));
                y=obj.lat_to_y(obj.lat(lat_id));

                overlay_data_grid(n,3)=obj.data(y,x);

           end
          
       end
       
       % plot zonal mean
       % varargin - NameValue Arguments of Line Properties
       function plot_zonal_mean(obj,varargin)
           
           zonal=zeros(obj.n_lat,1);
           for n=1:obj.n_lat
               zonal(n,1)=nanmean(obj.data(n,:));
           end
           
           %figure
           L=plot(obj.lat,zonal);
           xlabel('Latitude')
           ylabel([obj.longname,' (',obj.units,')'])
           box on;
           % set properties
           for n=1:2:numel(varargin)
               set(L,varargin{n},varargin{n+1});
           end
           

       end
       
       function [data,lon]=change_origin(obj)
           
           if isempty(obj.lon_origin)
               data=obj.data;
               lon=obj.lon;
           else
               new_lon_origin=obj.lon_origin;
               % expand longitudes out
               lon_left_edge=[-flipud(obj.lon+unique(diff(obj.lon))/2+360); obj.lon-unique(diff(obj.lon))/2; obj.lon-unique(diff(obj.lon))/2+360];
               tmp_data=[obj.data obj.data obj.data];
               
               % check new origin exists
                if ismember(new_lon_origin,lon_left_edge)
                    ind=find(lon_left_edge==new_lon_origin);
                else
                    data=obj.data;
                    lon=obj.lon;
                    disp('Longitude origin selected not available. Choose from:')
                    disp(lon_left_edge)
                    return
                end
                
                data=tmp_data(:,ind:ind+obj.n_lon-1);
                lon=lon_left_edge(ind:ind+obj.n_lon-1)+unique(diff(obj.lon))/2;
                
           end
           
               
           end
           
%            % check new origin exists
%            if ismember(new_lon_origin,lon_left_edge)
%                ind=find(lon_left_edge==new_lon_origin);
%            else
%                data=[]; lon=[];
%                disp('Longitude origin selected not available. Choose from:')
%                disp(lon_left_edge)
%                return
%            end
%            
%            
%            
%            if isempty(obj.lon_origin)
%                n=0;
%            elseif ismember(obj.lon_origin,lon_left_edge)
%                 n=find(lon_left_edge==obj.lon_origin);
%            else
%                data=[];
%                disp('Longitude origin selected not available. Choose from:')
%                disp(lon_left_edge)
%                return
%            end
%            
%            data=circshift(obj.data,n,2);
          
       
       % re-plot figure wwhen user changes properties
       function handlePropertyEvents(obj,src,evnt)
          if any(src.Name) & obj.autoplot
                evnt.AffectedObject.plot
          end
       end
       
        
    end
    
    methods (Static)
       
       % save png
       function save_bitmap(filename)
           set(gcf, 'InvertHardcopy', 'off')
           print(filename,'-dpng')
       end
       
       % save svg
       function save_vector(filename)
           print(filename,'-dsvg')
       end
       
       
       function [FCMAP] = make_cmap(PCNAME,PNCOL)
        % make_cmap
        %
        %   *******************************************************************   %
        %   *** make colorbar color scale *************************************   %
        %   *******************************************************************   %
        %
        %   make_cmap( ... )
        %   creates an evenly-spaced color scale of a given number of colors
        %
        %   PCNAME [STRING]
        %   --> name of colorbar colorscale
        %   PNCOL [integer]
        %   --> the number of colors
        %
        %   The color scale options are:
        %
        %   (1) One of the original MATLAB colorscales:
        %       'jet','hsv','hot','cool','spring','summer','autumn','winter',
        %       'gray','bone','copper','pink'
        %   (2) The new MATLAB (2014 release) default (the 'jet' replacement):
        %       'parula'
        %       (note that this color scheme is only approximated here)
        %   (3) An alternative rainbow replacements:
        %       'CubicYF'
        %       'LinearL'
        %       http://mycarta.wordpress.com/2013/02/21/perceptual-rainbow-palette-
        %       the-method/
        %       http://www.mathworks.com/matlabcentral/fileexchange/28982-perceptua
        %       lly-improved-colormaps/content/pmkmp/pmkmp.m
        %       (note that this color scheme is only approximated here)
        %   (4) One of the divergent cbrewer schemes
        %       'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'
        %       http://colorbrewer2.org/
        %       http://uk.mathworks.com/matlabcentral/fileexchange/34087-cbrewer---
        %       colorbrewer-schemes-for-matlab
        %       (note that these color schemes are only approximated here)
        %   (5) One of the built-in/internal color schemes:
        %       'anom' -- is the existing anomoly color scehem (blue-white-red)
        %       '_parula' -- is an approximation of the parula m-function
        %
        %   The default is 'parula', and if parula.m is not available,
        %   the approximateion '_parula' is used instead.
        %
        %   Original author:
        %   Andy Ridgwell <andy@seao2.org>
        %
        %   License:
        %   CC BY-NC-SA 3.0 (http://creativecommons.org/licenses/by-nc-sa/3.0/)
        %
        %   *******************************************************************   %

        % *********************************************************************** %
        % ***** HISTORY ********************************************************* %
        % *********************************************************************** %
        %
        %   14/12/03: CREATED
        %   15/01/01: adjustment to default colotmap
        %   15/01/09: re-wrote the 're-gridding' of the color map
        %   15/01/11: added 'LinearL' color scheme
        %   15/01/11: renamed from make_cmap5.m
        %             + substituted simple linear interpolation in place of 
        %              existing bizarre and pooor code ...
        %   15/01/11: added wt% color scale
        %   15/02/11: added missing scaling to anom and wt% color scales
        %   15/03/01: changed 'cubic' -> 'PCHIP' in interp1 on MATLAB 'advice' ...
        %   15/03/03: sorted out the previous 'fix' that was bugged on old versions
        %             improved warning message
        %   16/02/17: corrected error message
        %   17/05/15: added white-to-red-black, white-to-blue-black scale options
        %
        % *********************************************************************** %

        % *********************************************************************** %
        % *** INITIALIZE PARAMETERS & VARIABLES ********************************* %
        % *********************************************************************** %
        %
        % process dummy parameters
        str_name = PCNAME;
        col_n = PNCOL;
        %
        if (col_n < 2),
            disp(['ERROR: number of colors (here = ', num2str(col_n), ') must be > 1.']);
            return;
        end
        % initialize color map array
        loc_str_cmap = [];
        % determine MUTLAB version
        tmp_mutlab = version('-release');
        str_mutlab = tmp_mutlab(1:4);
        par_mutlab = str2num(str_mutlab);
        %
        % *********************************************************************** %

        % *********************************************************************** %
        % *** CREATE COLOR SCALE ************************************************ %
        % *********************************************************************** %
        %
        % *** SET DEFINING TIE-POINT COLORS ************************************* %
        %
        switch str_name
            case {'anom'}
                loc_cmap = [  0   0 255;
                            127 255 255;
                            255 255 255;
                            255 255 127;
                            255   0   0];
                loc_cmap = loc_cmap/255.0;  
            case {'wt%'}
                loc_cmap = [ 64  32  16;
                            128  64  32;
                            255 128  64;
                            255 255 255;
                            255 255 255];
                loc_cmap = loc_cmap/255.0;  
            case {'red'}
                loc_cmap = [255 255 255;
                            255 255 127;
                            255   0   0;
                              0   0   0];
                loc_cmap = loc_cmap/255.0;  
            case {'blue'}
                loc_cmap = [255 255 255;
                            127 255 255;
                              0   0 255;
                              0   0   0];
                loc_cmap = loc_cmap/255.0;  
            case {'_parula'}
                loc_cmap = [0.208100000000000   0.166300000000000   0.529200000000000;
                            0.079475000000000   0.515900000000000   0.832825000000000;
                            0.198550000000000   0.721400000000000   0.630950000000000;
                            0.826575000000000   0.732025000000000   0.346350000000000;
                            0.976300000000000   0.983100000000000   0.053800000000000];
            case {'RdYlGn'}
                loc_cmap = [                0   0.407843137254902   0.215686274509804;
                            0.525490196078432   0.796078431372549   0.401960784313726;
                            1.000000000000000   1.000000000000000   0.749019607843137;
                            0.974509803921569   0.554901960784314   0.321568627450980;
                            0.647058823529412                   0   0.149019607843137];
            case {'RdYlBu'}
                loc_cmap = [0.192156862745098   0.211764705882353   0.584313725490196;
                            0.562745098039216   0.764705882352942   0.866666666666666;
                            1.000000000000000   1.000000000000000   0.749019607843137;
                            0.974509803921569   0.554901960784314   0.321568627450980;
                            0.647058823529412                   0   0.149019607843137];
            case {'RdGy'}
                loc_cmap = [0.101960784313725   0.101960784313725   0.101960784313725;
                            0.629411764705882   0.629411764705882   0.629411764705882;
                                            1                   1                   1;
                            0.898039215686275   0.511764705882353   0.405882352941176;
                            0.403921568627451                   0   0.121568627450980];
            case {'RdBu'}
                loc_cmap = [0.019607843137255   0.188235294117647   0.380392156862745;
                            0.417647058823529   0.674509803921569   0.817647058823529;
                            0.968627450980392   0.968627450980392   0.968627450980392;
                            0.898039215686275   0.511764705882353   0.405882352941176;
                            0.403921568627451                   0   0.121568627450980];
            case {'PuOr'}
                loc_cmap = [0.176470588235294                   0   0.294117647058824;
                            0.600000000000000   0.560784313725490   0.749019607843137;
                            0.968627450980392   0.968627450980392   0.968627450980392;
                            0.935294117647059   0.615686274509804   0.233333333333333;
                            0.498039215686275   0.231372549019608   0.031372549019608];
            case {'PRGn'}
                loc_cmap = [                0   0.266666666666667   0.105882352941176;
                            0.501960784313725   0.770588235294118   0.503921568627451;
                            0.968627450980392   0.968627450980392   0.968627450980392;
                            0.680392156862745   0.543137254901961   0.741176470588235;
                            0.250980392156863                   0   0.294117647058824];
            case {'PiYG'}
                loc_cmap = [0.152941176470588   0.392156862745098   0.098039215686275;
                            0.609803921568628   0.809803921568627   0.390196078431373;
                            0.968627450980392   0.968627450980392   0.968627450980392;
                            0.907843137254902   0.590196078431373   0.768627450980392;
                            0.556862745098039   0.003921568627451   0.321568627450980];
            case {'BrBG'}
                loc_cmap = [                0   0.235294117647059   0.188235294117647;
                            0.354901960784314   0.698039215686274   0.658823529411765;
                            0.960784313725490   0.960784313725490   0.960784313725490;
                            0.811764705882353   0.633333333333333   0.333333333333333;
                            0.329411764705882   0.188235294117647   0.019607843137255];
            case {'LinearL'}
                loc_cmap =  [0.0143	0.0143	0.0143;
                             0.1413	0.0555	0.1256;
                             0.1761	0.0911	0.2782;
                             0.1710	0.1314	0.4540;
                             0.1074	0.2234	0.4984;
                             0.0686	0.3044	0.5068;
                             0.0008	0.3927	0.4267;
                             0.0000	0.4763	0.3464;
                             0.0000	0.5565	0.2469;
                             0.0000	0.6381	0.1638;
                             0.2167	0.6966	0.0000;
                             0.3898	0.7563	0.0000;
                             0.6912	0.7795	0.0000;
                             0.8548	0.8041	0.4555;
                             0.9712	0.8429	0.7287;
                             0.9692	0.9273	0.8961]; 
            case {'CubicYF'}
                loc_cmap =  [0.5151    0.0482    0.6697;
                             0.5199    0.1762    0.8083;
                             0.4884    0.2912    0.9234;
                             0.4297    0.3855    0.9921;
                             0.3893    0.4792    0.9775;
                             0.3337    0.5650    0.9056;
                             0.2795    0.6419    0.8287;
                             0.2210    0.7123    0.7258;
                             0.2468    0.7612    0.6248;
                             0.2833    0.8125    0.5069;
                             0.3198    0.8492    0.3956;
                             0.3602    0.8896    0.2919;
                             0.4568    0.9136    0.3018;
                             0.6033    0.9255    0.3295;
                             0.7066    0.9255    0.3414;
                             0.8000    0.9255    0.3529]; 
            case {'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink'}
                if (exist(str_name) == 0),
                    disp(['ERROR: colormap ', str_name, ' does not exist.']);
                    return;
                else
                    loc_str_cmap = str_name;
                end
            case {'parula'}
                if (exist('parula') == 2 || exist('parula') == 5),
                    loc_str_cmap = 'parula';
                else
                    loc_cmap = [0.208100000000000   0.166300000000000   0.529200000000000;
                                0.079475000000000   0.515900000000000   0.832825000000000;
                                0.198550000000000   0.721400000000000   0.630950000000000;
                                0.826575000000000   0.732025000000000   0.346350000000000;
                                0.976300000000000   0.983100000000000   0.053800000000000];
                    disp(['WARNING: colormap ', str_name, ' cannot be found (MATLAB version: ', num2str(par_mutlab), ' too old) => using parula-like replacement.']);
                end
            otherwise
                    loc_cmap = [0.208100000000000   0.166300000000000   0.529200000000000;
                                0.079475000000000   0.515900000000000   0.832825000000000;
                                0.198550000000000   0.721400000000000   0.630950000000000;
                                0.826575000000000   0.732025000000000   0.346350000000000;
                                0.976300000000000   0.983100000000000   0.053800000000000];
                disp(['WARNING: colormap ', str_name, ' cannot be found => using parula-like default.']);
        end
        %
        % *** CREATE COLOR SCALE ************************************************ %
        %
        if ~isempty(loc_str_cmap),
            cmap = colormap(eval([loc_str_cmap '(' num2str(col_n) ')']));
        else
            loc_n_old = [1:1:length(loc_cmap(:,1))];
            loc_n_new = [1:((length(loc_cmap(:,1))-1)/(col_n-1)):length(loc_cmap(:,1))];
            if (par_mutlab >= 2014),
                cmap = interp1(loc_n_old(:),loc_cmap(:,:),loc_n_new(:),'pchip');
            else
                cmap = interp1(loc_n_old(:),loc_cmap(:,:),loc_n_new(:),'cubic');
            end
        end
        %
        % *********************************************************************** %

        % *********************************************************************** %
        % *** END *************************************************************** %
        % *********************************************************************** %
        %
        % export data
        FCMAP = cmap;
        %
        % *********************************************************************** %
        
       end

       
       
    end
        


    
end