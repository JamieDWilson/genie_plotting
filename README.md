# genie_plotting
Quick interactive plots of lon/lat or lat/depth fields from genie output.

---

<u>Contents<u>
- [Inputs and Outputs](#inputs-and-outputs)
- [Basic Use](#basic-use)
- [Interactive Plotting](#interactive-plotting)
- [Static Plotting and Subplots](#static-plotting-and-subplots)
- [Plotting custom arrays](#plotting-custom-arrays)
- [Depth Integrated](#depth-integrated)
- [Zonal Averages](#zonal-averages)
- [Adding Contours](#adding-contours)
- [Overlaying datapoints](#overlaying-datapoints)
- [Saving plots](#saving-plots)
- [Misc Tips](#misc-tips)
- [Known Issues](#known-issues)
- [Features to add](#features-to-add)


---
  
## Inputs and Outputs
Files: 
- plot_genie_lonlat: plots longitude/latitude field from either 3D or 2D output. 
- plot_genie_latdepth: plots latitude/depth field from 3D output
- plot_genie: contains default parameter values and loads data

 Inputs: 
- `output_dir`: name or path of genie output file, e.g., 'fields_biogem_3d.nc' or 'experiment_dir/biogem/fields_biogem_3d.nc'
- `var_name`: name of variable in netCDF to plot, e.g., 'ocn_PO4'
- `year`: timeslice to plot, e.g., 9999.5
- `depth` or `longitude`: level to plot (for depth, 1 for surface, >1 for interior)
- `data_scale`: scale output by orders of magnitude, e.g., 1e-6 (micro)
- `optional #1`: name or path to 2nd output file to plot as difference
- `optional #2`: timeslice of 2nd output to plot as difference

Outputs:
- `object handle`: used to interactively change the plot and/or call subfunctions

---

## Basic Use

Plot surface PO$_4$ lon/lat field from 3D netcdf output and scale to micromolar units:
```matlab
fig=plot_genie_lonlat( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 );
```

Plot surface PO$_4$ lon/lat field as difference with another experiment. The colourscale will adjust to a blue-white-red difference scale:
```matlab
fig=plot_genie_lonlat( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 , 'fields_biogem_3d_alt.nc' , 9999.5 );
```

If a variable name, year, or depth/longitude is not available, the functions will list available options.

---

## Interactive Plotting

By default, calling plot_genie_lonat() or plot_genie_latdepth() will make a plot in a new window and creates a handle for that plot. The handle can then be used to change various parameters and the plot will automatically re-draw itself. The handle refers to its specific plot so multiple plots can be opened and edited independently. 

Plot surface PO$_4$ lon/lat field and adjust colourscale:
```matlab
fig=plot_genie_lonlat( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 );
fig.cmin=0;
fig.cmax=2;
fig.c_n_levels=10;
fig.title_text='Surface PO4';
```

The following parameters are editable and interactive:

*Colourscale and Colourbar*
- `cmin`: minimum value for colorscale (float, []=min of data)
- `cmax`: maximum value for colorscale (float, []=max of data)
- `c_nlevels`: number of colour intervals (integer)
- `colormap`: name of colormap (string)
- `reverse_colormap`: reverse the colormap? (logical)
- `add_colorbar`: plot a colourbar (logical)
- `colorbar_ticks`: set ticks for colorbar (float. []=auto)
- `colorbar_text`: colorbar text (string, 'auto' = uses variable name and units from netcdf file)

*Map Display*
- `lon_label`: have a longitude label? (logical)
- `lat_label`: have a latitude label? (logical)
- `z_label`: have a depth label? (logical)
- `lat_tick_label`: latitude tick labels to display in deg N (float)
- `lon_tick_label`: longitude tick labels to display in deg E(float)
- `z_tick_label`: depth tick labels to display in km (float)
- `lon_origin`: longitude of origin for plot in deg E (float, []=origin in netcdf file)
- `lat_limits`: set latitude limits in deg N (float)
- `lon_limits`: set longitude limits in deg E (float)

*Figure Display*
- `title_text`: text for plot title (string, 'auto'=automatic using filename, variable name, depth/longitude and year)

---

## Static Plotting and Subplots

You can turn off interactive plotting by setting `autoplot=false` in the plot_genie file. Note, the figure handle now points to the current open figure window!

Parameters can be changed and a final plot created by calling the .plot function. This is useful if you want a create a fully formatted plot in a script:
```matlab
fig=plot_genie_lonlat( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 );
fig.cmin=0;
fig.cmax=2;
fig.c_n_levels=10;
fig.plot % plot the final figure
```

This functionality also allows multi-panel figures to be constructed via the subplot command:
```matlab
figure
subplot(2,1,1)
fig=plot_genie_lonlat ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 4999.5 , 1 , 1e-6 );
fig.cmin=0;
fig.cmax=2;
fig.c_n_levels=10;
fig.plot;
subplot(2,1,2)
fig=plot_genie_lonlat( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 );
fig.cmin=0;
fig.cmax=2;
fig.c_n_levels=10;
fig.plot;
```
---

## Plotting Custom Arrays

You can plot an output derived from the netcdf file outputs by replacing the `var_name` input with a matlab array. The array needs to have NaNs for land and be orientated such that imagesc(array) looks the right way up/round. The plot command still requires a relevant netcdf file to read in the grid information and the `data_scale` parameter will still apply.

An example of the workflow:

```matlab
% read in raw netcdf data
data=ncread('fields_biogem_3d.nc','ocn_PO4');
% slice last year and surface
data=data(:,:,1,end);
% set missing values to nans
data(data>1e29)=NaN;
% orientate 
data=rot90(data);
imagesc(data) % quick visual check
% do something - e.g., percentage of total
data=data./nansum(reshape(data,36*36,1))*100;

% plot
fig=plot_genie_lonlat('fields_biogem_3d.nc',data,1,1e-6);
```

---

## Depth Integrated

Not yet implemented!

---

## Zonal Averages

plot_genie_latdepth can average zonally by setting a range of longitudes in the function call. This call creates a global zonal mean:
```matlab
fig=plot_genie_latdepth ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , [1:36] , 1e-6 );
```

plot_genie_lonlat can also create a separate zonal average plot and output the data:
```matlab
fig=plot_genie_lonlat ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 );
[ latitude , zonal_mean ]=fig.plot_zonal_mean;
```

The zonal average plot can be formatted as a standard matlab plot by passing [Line Properties](https://www.mathworks.com/help/matlab/ref/matlab.graphics.chart.primitive.line-properties.html):
```matlab
fig=plot_genie_lonlat ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 );
[ latitude , zonal_mean ] fig.plot_zonal_mean('LineWidth',1.0,'LineStyle',':');
```

---

## Adding Contours

Contours can be added to a plot using the .add_contours function. The takes the same arguments as the contour() function in MATLAB: levels is the number of contours or the specified levels, and name-value pair arguments can be passed for formatting:

```matlab
fig=plot_genie_lonlat ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 );
fig.cmin=0;
fig.cmax=2;
fig.add_contours([0:0.2:2],'Color','white');
```

Note, this currently works better with the lon/lat plots. The lat/depth plots produce very blocky contours!


---

## Overlaying datapoints

Datapoints can be overlain on the plots with colours matching the colourscale. The overlay data needs to be in an nx3 array of longitude, latitude, data:
```matlab
% create overlay data array
data_lon=[23; -180; -2.3; 80]; % deg E
data_lat=[-30; -22.2; 22.2; 30]; % deg N
data_val=rand(4,1); 
data_overlay=[data_lon data_lat data_val]; % into one array

fig=plot_genie_lonlat ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 );
fig.overlay_data(data_overlay);
```
Set the data to NaNs to plot a grey circle instead of coloured. This is useful to display locations of data. 

The overlay data point size can be altered via the interactive variable: `overlay_point_size`, e.g., `fig.overlay_point_size=60;'.

You can also extract the equivalent genie data from the nearest grid-cell to an output array with grid longitde, grid latitude, and grid data:
```matlab
fig=plot_genie_lonlat ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 );
[grid_values]=fig.overlay_data(data_overlay);
```

---

## Saving plots

Any plot can be saved to a bitmap (png) or vector format (svg) using the .save_bitmap or .save_vector functions via the figure handle:
```matlab
fig=plot_genie_lonlat ( 'fields_biogem_3d.nc' , 'ocn_PO4' , 9999.5 , 1 , 1e-6 );
fig.save_bitmap('genie_PO4.png')
fig.save_vector('genie_PO4.svg')
```

---

## Misc Tips

*Colourscale Options*
- Matlab in-built: 'jet', 'hsv', 'hot', 'cool', 'spring', 'summer', 'autumn', 'winter', 'gray', 'bone', 'copper', pink', 'parula'
- Rainbow replacements: 'CubicYF', 'LinearL'
- Divergent: 'BrBG', 'PiYG', 'PRGn', 'PuOr', 'RdBu', 'RdGy', 'RdYlBu', 'RdYlGn'
- Thanks to Andy Ridgwell!  

*Plotting with 2D netcdf variables*
- the `depth` input to plot_genie_lonlat is used only for making the title
- ecogem output seems to use the depth numbering backwards so setting `depth` to 16 rather than 1 may state the correct depth in the title
- variables like grid_topo in biogem output are sometimes referenced to the year 0.5

----

## Known Issues

- overlaid data may not be displayed if the lon/lat values are close to the edge of the figure edges.
- the overlay data function has to be called after fig.plot if not using the functions interactively.
- the figure background is set to the grey continent colour.
- contouring a lat/depth plot is very blocky (due to using the imagesc function for plotting)

---

## Features to add
  
- depth-integration
- circulation vectors
- reference lists for depth/longitude and colourscale
- custom relationships between two datasets



