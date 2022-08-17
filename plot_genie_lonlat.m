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



classdef plot_genie_lonlat < plot_genie
    
    methods
        
        function obj=plot_genie_lonlat ( output_file , var_name , year , depth , data_scale , varargin )
            
            
            obj=obj@plot_genie( output_file , var_name , year , depth , data_scale );
            obj.opt_args=varargin;
            obj.ind_var_name='depth';
            obj.depth=obj.ind_var;

            % load in plot data
            initialise_plot_data(obj);
            
            % initial plot
            if obj.autoplot
                plot(obj);
            end
            
        end
        
    end
    

    methods(Access = public)
        
        function [] = plot(obj)
            
            % change origin
            [plot_data,plot_lon]=obj.change_origin;
            % keep copy of plot_data
            obj.plot_data=plot_data;
            obj.lon_limits=[min(plot_lon) max(plot_lon)];
            
            % make coastlines
            obj.make_coastlines(plot_data);
            
            % make colormap
            %obj.make_colormap;
            
            % point at relevant figure handle
            if obj.autoplot
                figure(obj.fig);
            end

            % plot data, NaNs not plotted
            if obj.autoplot
                clf
            end
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
            obj.make_colormap;
            colormap(obj.ax,obj.c);
            caxis(obj.ax,[obj.cmin obj.cmax]);
            
            % colorbar
            obj.plot_colorbar;
            
            % title
            obj.auto_title_text;
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
            
            % apply lon/lat plot limits
            set(obj.ax,'XLim',obj.lon_to_x(obj.lon_limits)+[-0.5 0.5],'YLim',fliplr(obj.lat_to_y(obj.lat_limits))+[-0.5 0.5]);
            
            
            
        end
            
        function [] = lonlat_to_xy(obj,lon,lat)
            
            obj.lon_to_x=griddedInterpolant(lon,[1:1:obj.n_lon]);
            obj.lat_to_y=griddedInterpolant(lat,[obj.n_lat:-1:1]);
            
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
       function [lat , zonal] = plot_zonal_mean(obj,varargin)
           if obj.autoplot
               figure;
           end
           
           % mean across latitudes, flipped to match latitude order
           zonal=flipud(nanmean(obj.data,2));
           
           %figure
           L=plot(obj.lat,zonal);
           xlabel('Latitude')
           ylabel([obj.longname,' (',obj.units,')'])
           box on;
           % set properties
           for n=1:2:numel(varargin)
               set(L,varargin{n},varargin{n+1});
           end
           lat=obj.lat;
           

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
                    disp(obj.lon-unique(diff(obj.lon))/2)
                    return
                end
                
                data=tmp_data(:,ind:ind+obj.n_lon-1);
                lon=lon_left_edge(ind:ind+obj.n_lon-1)+unique(diff(obj.lon))/2;
                
           end
           
               
       end
       
       function [] = auto_title_text(obj)
           
           % auto title text
            if strmatch(obj.title_text','auto')
                if obj.data_2_flag
                    obj.title_text={ [obj.longname ' at ' num2str(round(obj.z(obj.depth),2)) ' m'] , ['(' obj.output_dirs{1} ' @ ' num2str(obj.time{1}(obj.year(1))) ' years) - (' obj.output_dirs{2} ' @ ' num2str(obj.time{2}(obj.year(2))) ' years)'] };
                else
                    obj.title_text={ [obj.longname ' at ' num2str(round(obj.z(obj.depth),2)) ' m'] , ['(' obj.output_dirs{1} ' @ ' num2str(obj.time{1}(obj.year(1))) ' years)'] };
                end
            end

           
       end
     
    end
    
        


    
end