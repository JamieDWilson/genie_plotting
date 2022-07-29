%% plot_genie_latdepth
% -------------------------------------------------------------------------
% Description: Quick interactive plots of lat/depth fields from genie output.
%
% Inputs: 
%   - output_dir: name or path of genie output file, e.g., 'fields_biogem_3d.nc' or 'experiment_dir/biogem/fields_biogem_3d.nc'
%   - var_name: name of variable in netCDF to plot, e.g., 'ocn_PO4'
%   - year: timeslice to plot, e.g., 9999.5
%   - latitude: latitude number to plot, e.g., 1 to 36. Can be specificed as a range, e.g., [1:36], and the plot will show the average
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
%   -> plot basic PO4 lat/depth field and adjust colourscale interactively
%   fig=plot_genie_latdepth ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 )
%   fig.cmin=0;
%   fig.cmax=2;
%   fig.c_n_levels=10;
%
%   -> plot difference with another experiment
%   fig=plot_genie_latdepth ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 , 'fields_biogem_3d_alt.nc' , 9999.5 )
%
%   -> plot zonal average PO4 lat/depth field and adjust colourscale interactively
%   fig=plot_genie_latdepth ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , [1:36] , 1e-6 )
%   fig.cmin=0;
%   fig.cmax=2;
%   fig.c_n_levels=10;
%
%   --------------------- STATIC PLOTTNG & SUBPLOTs -----------------------
%
%   - setting autoplot=false allows different usage. Parameter values can
%   be set and the plot drawn finally using the .plot function. This allows
%   plots to be drawn as subplots. This draws to the currently
%   selected figure window.
%
%   -> plot single final figure, e.g., in a script
%   fig=plot_genie_latdepth ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 )
%   fig.cmin=0;
%   fig.cmax=2;
%   fig.c_n_levels=10;
%   fig.plot;
%   
%   -> plot as a subplot
%   figure
%   subplot(2,1,1)
%   fig=plot_genie_latdepth ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 )
%   fig.cmin=0;
%   fig.cmax=2;
%   fig.c_n_levels=10;
%   fig.plot;
%   subplot(2,1,2)
%   fig=plot_genie_latdepth ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 )
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
%   - no overlay data function yet
% -------------------------------------------------------------------------

classdef plot_genie_latdepth < plot_genie
    
    methods
        
        function obj = plot_genie_latdepth ( output_file , var_name , year , longitude , data_scale , varargin )
            
            obj=obj@plot_genie( output_file , var_name , year , longitude , data_scale );
            obj.opt_args=varargin;
            obj.ind_var_name='lon';
            obj.longitude=obj.ind_var;
            
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
            
            % expand out depth axis for imagesc
            [plot_data,plot_z]=obj.expand_z;
            % keep copy of plot_data
            obj.plot_data=plot_data;
            
            % make coastlines
            obj.make_coastlines(plot_data);
            
            % make colormap
            obj.make_colormap;
            if obj.reverse_colormap
                obj.c=flipud(obj.c);
            end   
            
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
            obj.im=imagesc(plot_data,'AlphaData',imAlpha);
            
            % get axes
            if ~obj.autoplot
                obj.ax=gca;
            else
                obj.ax=obj.fig.CurrentAxes;
            end
            
            % lon/lat axes formatting
            obj.latdepth_to_xy(obj.lat,plot_z);
            obj.plot_ticks;
            set(obj.ax,'YTick',obj.z_ticks,'YTickLabel',num2cell(obj.z_tick_label),'XTick',obj.lat_ticks,'XTickLabel',num2cell(obj.lat_tick_label)) % set custom ticks/labels
            if obj.lat_label
                xlabel('Latitude')
            end
            if obj.z_label
                ylabel('Depth (km)')
            end
            
            % colorbar
            if obj.colorbar
                C=colorbar;
                ylabel(C,obj.colorbar_text);
            end
            
            % plot formating
            axis normal
            pbaspect([1.5 1 1]) % rectangular aspect
            hold on
    
            % plot coastlines
            for i=1:numel(obj.coastlines(:,1))
                plot(obj.coastlines(i,1:2)-1.5,obj.coastlines(i,3:4)-1.5,'k')
                hold on
            end
            
            % set background (and land) color
            set(obj.ax, 'color', [0.8 0.8 0.8]);
            
            % colormaps
            colormap(obj.ax,obj.c);
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
        
        
        function [datai,zi_ind] = expand_z(obj)
            
            % expand data array to approximate correct grid sizes
            % find nearest integer ratio of each depth level to surface, expand grid
            % by ratio - is approximate but I hate pcolor.
            dz=diff(obj.zt_edges);
            dz_ratio_surface=round(dz./dz(1));
            n_zi=sum(dz_ratio_surface);
            datai=zeros(n_zi,obj.n_lat);
            zi_ind=zeros(obj.n_z,1);

            pos=1;
            for n=1:obj.n_z
                for nn=1:obj.n_lat
                    ind=[pos:pos+dz_ratio_surface(n)-1];
                    datai(ind,nn)=obj.data(n,nn);
                end
                zi_ind(n,1)=(min(ind)+max(ind))/2;
                pos=pos+numel(ind);

            end
            
            
        end
        
        
        function [] = latdepth_to_xy(obj,lat,z)
            
            obj.lat_to_x=griddedInterpolant(obj.lat,[1:1:obj.n_lat]);
            obj.z_to_y=griddedInterpolant(obj.z,z);
            
        end
        

        
        function [] = plot_ticks(obj)
            
            obj.lat_ticks=obj.lat_to_x(obj.lat_tick_label);
            obj.z_ticks=obj.z_to_y(obj.z_tick_label*1000);

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
        
    end
    
end