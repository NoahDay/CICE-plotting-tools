function [h,cb] = bedmap2(mapvar,varargin)
% bedmap2 uses simple commands to plot Bedmap2 data. If you
% have Matlab's Mapping Toolbox, data are plotted in georeferenced
% coordinates. If you do not have the Mapping Toolbox, data are
% automatically plotted in polar stereographic coordinates. 
% 
%% Requirements
% This function requires The Antarctic Mapping Tools package available 
% for free on the Mathworks File exchange website here: 
% http://www.mathworks.com/matlabcentral/fileexchange/47638-antarctic-mapping-tools
% 
%% Syntax
% 
%  bedmap2('mapvar') 
%  bedmap2(...,'PropertyName',PropertyValue,...) 
%  bedmap2(...,'res',ResolutionKilometers) 
%  bedmap2(...,'xy') 
%  bedmap2(...,'xy','km') 
%  bedmap2(...,'colorbar','ColorbarLocation') 
%  bedmap2(...,'frame','on') 
%  bedmap2(...,'rotation',RotationDegrees) 
%  bedmap2(...,'oceancolor',ColorSpec) 
%  [h,cb] = bedmap2(...)
% 
%% Description 
% 
% bedmap2('mapvar') plots a continent-wide map of 'mapvar', which can
% be any of the following Bedmap2 variables: 
% 
% * 'gl'             outlines grounding line with optional linestyle arguments.
% * 'patchgl'        creates patches for the 200 largest grounded regions.
% * 'shelves'        outlines ice shelves with optional linestyle arguments. 
% * 'patchshelves'   creates patches for the 200 largest floating ice shelves.
% * 'coast'          outlines Antarctica's coast line with optional linestyle arguments. 
% * 'patchcoast'     creates patches for the 70 largest coast-bound regions.
% * 'thickness'      plots ice thickness with pcolorm, (units: m). 
% * 'bed'            plots bed elev. with pcolorm, (m relative to gl04c geoid). 
% * 'bedw            plots bed elev. with pcolorm, (m relative to wgs84 ellipsoid). 
% * 'bedc'           contour map of bed elevations given by zval (gl04c geoid).
% * 'bedcw'          contour map of bed elevations given by zval (wgs84 ellipsoid).
% * 'surface'        plots surface elev. with pcolorm, (m relative to gl04c geoid). 
% * 'surfacew        plots surface elev. with pcolorm, (m rel. to wgs84 ellipsoid). 
% * 'surfc'          contour map of surface elevations given by zval (gl04c geoid).
% * 'surfcw'         contour map of surface elevs. given by zval (wgs84 ellipsoid).
% * 'thickc'         contour map of ice thickness given by zvals. 
% * 'beduncertainty' plots the given uncertainty of bed elevation, (units: m).
% * 'geoid'          plots surface of Gl04c geoid (z_wgs84 = z_bedmap2 + gl04c)   
% * 'coverage'       binary grid showing distribution of ice thickness data used in the grid of ice thickness 
% * 'icemask'        0 = grounded, 1 = ice shelf, 127 = ocean
% * 'rockmask'       binary indicating presence of rocks
% * 'vostok'         binary indicating presence of Lake Vostok
% * ''               Calling bedmap2 with an empty string initializes a map via antmap. 
% 
% bedmap2(...,'PropertyName',PropertyValue,...) formats line or patch
% objects (e.g., 'LineWidth','Facealpha',etc.
% 
% bedmap2(...,'res',ResolutionKilometers) specifies a resolution for
% gridded raster data set plotting. By default, resolution is set
% dynamically using extents of current map if a map is already open. For
% example, a continent-wide map of bed elevation is downsampled to 5 km
% by default, but if a map is inialized and zoomed to the extents of a
% glacier, data may be plotted at 2 km or 1 km resolution. Plotting a
% continent-wide map of the bed at 1 km takes nearly a full minute for my
% computer to render, whereas the same geographic extents plotted at 5 km
% is rendered in 1.1 second. 
% 
% bedmap2(...,'xy') plots data in polar stereographic xy coordinates with
% units of meters (true latitude 71 S). By default, if you have a Mapping
% Toolbox license, data are plotted in georeferenced coordinates. If you do not
% have the Mapping Toolbox, the 'xy' option is automatically selected. 
% 
% bedmap2(...,'xy','km') plots in polar stereographic kilometers instead of 
% the default meters.  
%
% bedmap2(...,'colorbar','ColorbarLocation') specifies a colorbar
% location. For example, 'colorbar','north' places a colorbar at the top 
% of a map. Choose 'none' to plot data without a colorbar. 
%
% bedmap2(...,'frame','on') places a frame around the map. 
%
% bedmap2(...,'rotation',RotationDegrees) rotates view by
% RotationDegrees from standard projection. Folks who study the Ross Ice
% Shelf do this often, to put south at the bottom of the map. 
%
% bedmap2(...,'oceancolor',ColorSpec) specifies background color of axes; effectively, 
% the color of the ocean when plotting a patch of the continent. 
%
% [h,cb] = bedmap2(...) returns handle(s) h of plotted object(s) and
% colorbar handle cb. If plotted fields are contours, cb is the contour matrix. 
% 
%% Example 1: Outlines  
% Map a simple grounding line: 
% 
% bedmap2 gl
% 
% Add a fat magenta coast line: 
% 
% bedmap2('coast','color','m','linewidth',2) 
% 
%% Example 2: Patches 
% We'll overlay the following patches on a graticule called by antmap. 
% Plot a green patch around a fat blue grounding line, specify ice shelf
% face color with an RGB triplet, and make the ocean red:
% 
% antmap('graticule') 
% bedmap2('patchgl','facecolor','g','edgecolor','blue','linewidth',2)
% bedmap2('patchshelves','facecolor',[1 .8 .3],'oceancolor','r')   
% 
%% Example 3: Surfaces and contours
% Here we plot a surface elevation map. Then overlay thick red elevation
% contours at 2 km and 3 km elevation, and thin black lines every 200 m in between: 
% 
% bedmap2 surface 
% bedmap2('surfc',[2000 3000],'color','r','linewidth',2)
% bedmap2('surfc',2200:200:2800)
% 
%% Example 4: Combining AMT tools
% If a map is initialized before calling bedmap2, the bedmap2 function loads and 
% plots only data within the range of the map extents. Resolution is automatically chosen to
% plot somewhat quickly.
% 
% In this example we initialize a 1200 kilometer wide map of the Antarctic
% Peninsula, plot bed elevations as a surface, apply relief shading with
% shadem, then overlay semitransparent blue ice shelves: 
% 
% mapzoom('antarctic peninsula',1200,'inset','northeast') 
% bedmap2('bed','colorbar','none') 
% shadem(-18)
% bedmap2('patchshelves','facecolor','blue',...
%     'edgecolor','blue','facealpha',.5)
% 
%% Example 5: No Mapping Toolbox?  No Problem!
% If you do not have Matlab's Mapping Toolbox, the bedmap2 function will automatically 
% plot in polar stereographic coordinates. If you do have the Mapping
% Toolbox and you would like to plot in polar stereographic coordinates,
% include the 'xy' tag like this: 
% 
% bedmap2('gl','xy') 
% bedmap2('coast','xy','color','m','linewidth',2)
% xlabel('easting (m)') 
% ylabel('northing (m)')
% 
%% Changes in Version 4.0 (February 2015)
% This plotting function was rewritten for Bedmap2 version 4.0. 
% 
% The up side of the changes: 
%  * Plotting small regions at high resolution is much faster now,
%  * Matlab's Mapping Toolbox is no longer required for simple plotting, 
%  * More of the Bedmap2 data sets (coverage, masks) can now be plotted.
%
% The down side of the changes: 
% This function was starting to suffer from some feature creep, so I
% stopped supporting the following commands: 
% 
%  * grid         (now use antmap('graticule') instead)
%  * graticule    (now use antmap('graticule') instead)
%  * velocity     (use measures function instead)
%  * accumulation (no replacement)
%  * caxis        (follow bedmap2 with a call of caxis)
%  * colormap     (follow bedmap2 with a call of colormap)
%  * latlim       (zoom before calling bedmap2) 
%  * lonlim       (zoom before calling bedmap2)
% 
% How you'll have to change your ways: 
% For speed, the bedmap2 function now loads only enough data to fill the
% extents of any map that is current. If no map is open when you call
% bedmap2, or if a map is current but covers a large geographical area, 
% a downsampled Bedmap2 data set will be plotted. If you work in a specific
% region, use mapzoom before calling bedmap2.
% 
% VERSION 4.2 (September 2015): Changed contouring feature to accept any of Matlab's
% contouring syntax.  
% 
%% References
% If this function is useful for you, please cite the following: 
% 
% Fretwell, P., et al. "Bedmap2: improved ice bed, surface and thickness datasets for Antarctica." The Cryosphere 7.1 (2013).
% <http://dx.doi.org/10.5194/tc-7-375-2013> 
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% Computers & Geosciences. 104 (2017) pp.151-157. http://dx.doi.org/10.1016/j.cageo.2016.08.003
% 
%% Author Info
% This function and supporting documentation were written by Chad A. Greene
% of the University of Texas at Austin's Institute for Geophysics (UTIG). 
% Chad's website is http://www.chadagreene.com. You can email Chad at
% chad@chadagreene.com.  Handwritten correspondence is nice, too. 
% 
% This function is part of the Bedmap2 Toolbox for Matlab Version 4.0,
% written February 2015.     
% 
% See also bedmap2_download, bedmap2_data, mapzoom, antmap, scarlabel, and bedmap2_interp. 

%% Plot patches if no input arguments: 

if nargin==0
    bedmap2 'patchshelves'
    bedmap2 'patchgl'
    return
end

%% Input checks:

% assert(nargin>0,'The bedmap2 function requires at least one input. What variable are you trying to map?') 
assert(isnumeric(mapvar)==0,'The bedmap2 function wants a string as the first input. It looks like you''ve tried to enter something numeric.') 
assert(exist('bedmap2_bed.tif','file')==2,'Cannot find the bedmap2_bed.tif file. This might mean it''s in a folder where Matlab cannot find it. Check out the bedmap2_download script.') 

%% Set defaults: 

% If user has Mapping Toolbox, make a georeferenced map:  
if license('test','map_toolbox') 
    MakeGeoMap = true; 
else
    MakeGeoMap = false; 
end

% Colorbar location:
cbloc = 'EastOutside'; 
IncludeColorbar = true; 

rotation = 0; 
oceancolor = 'none'; 
frame = false; 


%% Parse inputs: 

if nargin>1
    
    % Even if user has a Mapping Toolbox license, they may want xy coordinates:  
    tmp = strcmpi(varargin,'xy'); 
    if any(tmp) 
        MakeGeoMap = false; 
        varargin = varargin(~tmp); 
    end
    

    % Check for user-defined colorbar location: 
    tmp = strncmpi(varargin,'cb',2)+strcmpi(varargin,'colorbar'); 
    if any(tmp) 
        cbloc = varargin{find(tmp)+1}; 
        tmp(find(tmp)+1) = 1; 
        varargin = varargin(~tmp); 
        if any(strcmpi({'off','none'},cbloc))
            IncludeColorbar = false; 
        end
    end
    
    % Contour elevation values? 
    tmp = strncmpi(varargin,'zlim',1)+strncmpi(varargin,'elev',4)+strcmpi(varargin,'isobath');
    if any(tmp)
        
        % Support the old style of defining z values: 
        if isscalar( varargin{find(tmp)+1})
            varargin{find(tmp)+1}= [varargin{find(tmp)+1} varargin{find(tmp)+1}];
        end
        varargin = varargin(~tmp); 
    end
        
    % Frame or box setting: 
    tmp = strcmpi(varargin,'box')+strcmpi(varargin,'frame'); 
    if any(tmp) 
        FrameSetting = varargin{find(tmp)+1}; 
        if any(strcmpi(FrameSetting,'on')) 
            frame = true; 
        end
        tmp(find(tmp)+1) = 1; 
        varargin = varargin(~tmp); 
    end
    
    % Ocean color:  
    tmp = strncmpi(varargin,'ocean',5);
    if any(tmp)
        oceancolor = varargin{find(tmp)+1}; 
        tmp(find(tmp)+1) = 1; 
        varargin = varargin(~tmp); 
        
        if MakeGeoMap % Because with a geo map, a frame must be there to give the ocean color.  
            frame = true; 
        end
    end
    
    % Rotation:  
    tmp = strncmpi(varargin,'rotation',3);
    if any(tmp)
        rotation = varargin{find(tmp)+1}; 
        tmp(find(tmp)+1) = 1; 
        varargin = varargin(~tmp); 
    end
    
end

% How many inputs are left? 
lv = length(varargin); 


%% Initialize a map or check limits of current map:  


% Get initial axis limits: 
if MakeGeoMap
    antmap % initializes map if necessary 
    
    xl = xlim; % get axes limits--should be roughly in the range of -0.48 to 0.48
    yl = ylim; 
    [lal,lol] = minvtran(xl,yl); % convert to lat,lon limits
    [xl,yl] = ll2ps(lal,lol); % convert to polar stereographic meters; 

else

    % If no axes already exist checking axis limits will initialize [0 1 0 1] axes. If this is the case, 
    % use the whole range of Bedmap2 data. Otherwise, if axes are already initialized, use current limits: 
    if sum(axis==[0 1 0 1])==4
        xl = 3333000*[-1 1];
        yl = 3333000*[1 -1]; 
    else
        xl = xlim; 
        yl = ylim; 
    end
end

%% Define data resolution: 

% Total area to plot in square kilometers: 
sqkm = (max(xl)-min(xl))*(max(yl)-min(yl)); 

% Set default resolution based on how much area will be plotted: 
res = interp1([1 3e12 5e12 1e13 2e13 9e15],[1 1 2 3 5 5],sqkm,'nearest');
        
% Overwrite automatic resolution calculation if user-defined: 
if lv>0
    tmp = strncmpi(varargin,'res',3); 
    if any(tmp)
        res = varargin{find(tmp)+1}; 
        tmp(find(tmp)+1) = 1; 
        varargin = varargin(~tmp); 
    end
end

%%  Load and plot data: 

switch lower(mapvar)
         
    case {'gl','grounding line','grounding','ground','gl outline',...
            'gloutline','outlinegl','outline gl'} 
        
        % Load data: 
        load('bedmap2gl.mat');
        
        % If user has not requested a line color, make it black: 
        if ~any(strcmpi(varargin,'color')+strncmpi(varargin,'linecolor',5)+strncmpi(varargin,'edgecolor',6))
            varargin{lv+1} = 'color'; 
            varargin{lv+2} = 'k'; 
        end
            
        % Plot data: 
        if MakeGeoMap
            [lat,lon] = polyjoin(gllat,gllon); 
            h = plotm(lat,lon,varargin{:}); 
        else
            h = plotps(cell2nancat(gllat),cell2nancat(gllon),varargin{:});
        end
        
        
    case {'coast','coast line','coastline','cl'}
        
        % Load data: 
        load('bedmap2coast.mat');
        
        % If user has not requested a line color, make it black: 
        if ~any(strcmpi(varargin,'color')+strncmpi(varargin,'linecolor',5)+strncmpi(varargin,'edgecolor',6))
            varargin{lv+1} = 'color'; 
            varargin{lv+2} = 'k'; 
        end
            
        % Plot coast line data: 
        if MakeGeoMap
            h = plotm(cell2nancat(coastlat),cell2nancat(coastlon),varargin{:});
        else
            h = plotps(cell2nancat(coastlat),cell2nancat(coastlon),varargin{:});
        end
        
        
    case {'patchgl','patch gl','fill gl','fillgl','gl patch','glpatch'}
        
        % Load data: 
        load('bedmap2gl.mat');
        
        % If user has not requested a face color, make it medium-ish gray: 
        if ~any(strcmpi(varargin,'color')+strcmpi(varargin,'facecolor'))
            varargin{lv+1} = 'facecolor'; 
            varargin{lv+2} = .6*[1 1 1];
            lv = lv+2; % Because we just lengthened varargin by 2. 
        end
        
        % If user has not requested an edge color, make it darkish gray: 
        if ~any(strcmpi(varargin,'edgecolor')+strcmpi(varargin,'linecolor'))
            varargin{lv+1} = 'edgecolor'; 
            varargin{lv+2} = .3*[1 1 1]; 
        end
        
        for k = 1:200 % plots the 200 largest grounded areas
            if MakeGeoMap
                h = patchm(gllat{k},gllon{k},varargin{:});
            else
                h = patchps(gllat{k},gllon{k},.6*[1 1 1],varargin{:});
            end
        end
        
        
    case {'patchcoast','patch coast','fill coast','fillcoast','coast patch',...
            'patchcoastline','patch coastline','patch coast line','coastpatch'}
        
        load('bedmap2coast.mat');
        
        % If user has not requested a face color, make it medium-ish gray: 
        if ~any(strcmpi(varargin,'color')+strcmpi(varargin,'facecolor'))
            varargin{lv+1} = 'facecolor'; 
            varargin{lv+2} = .6*[1 1 1];
            lv = lv+2; % Because we just lengthened varargin by 2. 
        end
        
        % If user has not requested an edge color, make it darkish gray: 
        if ~any(strcmpi(varargin,'edgecolor')+strcmpi(varargin,'linecolor'))
            varargin{lv+1} = 'edgecolor'; 
            varargin{lv+2} = .3*[1 1 1]; 
        end
        
        for k = 1:70 % plots the 70 largest coast-bound areas
            if MakeGeoMap
                h = patchm(coastlat{k},coastlon{k},varargin{:});
            else
                h = patchps(coastlat{k},coastlon{k},.6*[1 1 1],varargin{:});
            end
        end
        
        
    case {'patchshelves','patch shelves','fillshelves','fill shelves',...
            'shelf patch','shelfpatch','its okay to be shelfish'}
        
        load('bedmap2shelf.mat');   
        
        % If user has not requested a face color, make it light gray: 
        if ~any(strcmpi(varargin,'color')+strcmpi(varargin,'facecolor'))
            varargin{lv+1} = 'facecolor'; 
            varargin{lv+2} = .8*[1 1 1];
            lv = lv+2; % Because we just lengthened varargin by 2. 
        end
        
        % If user has not requested an edge color, make it dark gray: 
        if ~any(strcmpi(varargin,'edgecolor')+strcmpi(varargin,'linecolor'))
            varargin{lv+1} = 'edgecolor'; 
            varargin{lv+2} = .2*[1 1 1]; 
        end
        
        for k = 1:200 % plots the 200 largest ice shelves
            if MakeGeoMap
                h = patchm(shlflat{k},shlflon{k},varargin{:});
            else
                h = patchps(shlflat{k},shlflon{k},.8*[1 1 1],varargin{:});
            end
        end  
        
        
    case {'shelves','outline shelves','plotshelves','shelf outline',...
            'shelfoutline','outlineshelves'} 
        
        load('bedmap2shelf.mat');   
        
        % If user has not requested a line color, make it black: 
        if ~any(strcmpi(varargin,'color')+strncmpi(varargin,'linecolor',5)+strncmpi(varargin,'edgecolor',6))
            varargin{lv+1} = 'color'; 
            varargin{lv+2} = 'k'; 
        end
            
        % Plot data: 
        for k = 1:604
            if MakeGeoMap
                h = plotm(shlflat{k},shlflon{k},varargin{:});
            else
                h = plotps(shlflat{k},shlflon{k},varargin{:});
            end
        end
        
    case {'bedcontour','bed contour','bedc','bedwcontour','bedw contour','bedwc','bedcw','bedcontourw',...
                'bedcwgs','bedwgsc','bed wigs, see','surfacecontour','surface contour','surfc','sfzc','sfzcontour',...
             'sfz contour','sfc contour','sfccontour','surfacecontourw','surface contour w','surfcw','sfzcw','sfzcontourw',...
             'sfz contourw','sfc contourw','sfccontourw','surfacewcontour',...
             'surface w contour','surface contour wgs','surfwc','sfzwc',...
             'sfzwcontour','wgs sfzc','wgs surfc','surfc wgs','surfcwgs'...
             'sfz w contour','sfc w contour','sfcwcontour','thicknesscontour','thickness contour','thickc','thckc','thckcontour',...
             'thck contour'} 
        
        
        % Figure out exactly which data to load: 
        switch lower(mapvar) 
            case {'bedcontour','bed contour','bedc'}; 
                mv = 'bed'; 
            case {'bedwcontour','bedw contour','bedwc','bedcw','bedcontourw',...
                'bedcwgs','bedwgsc','bed wigs, see'}; 
                mv = 'bedw'; 
            case {'surfacecontour','surface contour','surfc','sfzc','sfzcontour',...
                    'sfz contour','sfc contour','sfccontour'}; 
                mv = 'surface'; 
            case {'surfacecontourw','surface contour w','surfcw','sfzcw','sfzcontourw',...
                    'sfz contourw','sfc contourw','sfccontourw','surfacewcontour',...
                    'surface w contour','surface contour wgs','surfwc','sfzwc',...
                    'sfzwcontour','wgs sfzc','wgs surfc','surfc wgs','surfcwgs'...
                    'sfz w contour','sfc w contour','sfcwcontour'}; 
                mv = 'surfacew'; 
            case {'thicknesscontour','thickness contour','thickc','thckc','thckcontour',...
                    'thck contour'}; 
                mv = 'thickness'; 
        end
                
        % Load data: 
        
        hold on;
            if MakeGeoMap
                [LAT,LON,Z] = bedmap2_data(mv,xl,yl,2*res,'resolution',res); 
                
                % Set NaNs to a near-zero finite value to ensure an existent coastline if user requests a zero-elevation contour:  
                if ~strncmpi(mapvar,'b',1)
                    Z(isnan(Z))=-0.01; 
                end
                
                [X,Y] = mfwdtran(LAT,LON); 
                [cb,h] = contour(X,Y,Z,varargin{:}); 
            else
                [X,Y,Z] = bedmap2_data(mv,xl,yl,2*res,'resolution',res,'xy'); 
                
                % Set NaNs to a near-zero finite value to ensure an existent coastline if user requests a zero-elevation contour:  
                if ~strncmpi(mapvar,'b',1)
                    Z(isnan(Z))=-0.01; 
                end
                [cb,h] = contour(X,Y,Z,varargin{:}); 
            end
        
        if ~MakeGeoMap
         daspect([1 1 1]) 
        end
                
        
    case '' 
        % Assume user just wanted to initialize a map; do nothing. 
        
    otherwise
        
        if MakeGeoMap
            [lat,lon,Z] = bedmap2_data(mapvar,xl,yl,2*res,'resolution',res); 
            
            if strcmpi(mapvar,'icemask')
                Z(Z==127) = 2; % allows 3 category classification.  
            end
            
            h = pcolorm(lat,lon,Z); 
            switch mapvar
                case {'bed','bottom','bedrock','depth','bedw','bottomw','bedrockw','depthw','bedwgs','bedwigs','bedwgs84',...
            'bed wgs','bed wgs 84'}
                    demcmap(Z,256);
                    colorbarstring = 'bed elevation (m)'; 
                    
                case {'surface','sfc','surf','sfz','surfacew','sfcw','surfw','sfzw','surfacewgs','sfcwgs','surfwgs','sfzwgs',...
                    'surfacewgs84','sfcwgs84','surfwgs84','sfzwgs84'}
                     demcmap(Z,256); 
                     colorbarstring = 'surface elevation (m)'; 
                    
                case {'beduncertainty','bed uncertainty','bedunc','beduncadunk'}
                    colormap(jet(16)) % Matches Fretwell 2013
                    caxis([66 700])   % Somewhat matches Fretwell 2013, but Fretwell et al. did not use a linear scale.  
                    colorbarstring = 'bed uncertainty (m)'; 
                    
                case  {'thickness','thick','th','thck'}
                    colorbarstring = 'ice thickness (m)'; 

                case {'coverage','cov'}
                    IncludeColorbar = false; % coverage is binary
                    colormap([1 0 0;1 1 1]); % white to red

                case {'icemask','ice'}       % has 3 categories: ocean, ice shelf, grounded.  
                    IncludeColorbar = false; 
                    colormap([0.5373    0.7412    0.4706;0.6039    0.4118    0.2039;1 1 1]); % three colors, ocean is white 
                    
                case {'rockmask','rock'} 
                    IncludeColorbar = false; 
                    colormap([1 0 0;1 1 1]); % white to red

                case {'lakemask','lake','vostok','vostokmask'} 
                    IncludeColorbar = false; 
                    colormap([1 0 0;1 1 1]); % white to red

                case {'geoid','gl04c','wgs','wgs84'}
                    colorbarstring = 'GL04c geoid relative to WGS84 ellipsoid (m)'; 
                    
                otherwise 
                    colorbarstring = []; 
            end
            
        else
            [x,y,Z] = bedmap2_data(mapvar,xl,yl,2*res,'resolution',res,'xy'); 
            h = pcolor(x,y,Z); 
            shading flat
            hold on
            colorbarstring = 'elevation (m)'; 
        end

        if IncludeColorbar
            cb = colorbar('location',cbloc); 
            switch lower(cbloc)
                case {'north','south','northoutside','southoutside'}
                    xlabel(cb,colorbarstring);
                otherwise
                    ylabel(cb,colorbarstring); 
            end
        end
        
        if ~MakeGeoMap
            daspect([1 1 1]) % Doing this after the colorbar placement tightens things up a bit. 
        end
                
end

%% Tinker with appearance: 

% Set rotation: 
view([rotation 90]); 

% Remove frame if requested: 
if MakeGeoMap
    if ~frame 
        set(gca,'visible','off'); % Removes the box frame if desired. 
        set(findall(gca, 'type', 'text'), 'visible', 'on');
    end
else
    box off
end

% Set ocean color:
set(gca,'color',oceancolor); 

%% Clean up: 

if nargout==0
    clear h; 
end


end



function B = cell2nancat(A) 
%cell2nancat concatenates elements of a cell into a NaN-separated vector. 
% 
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas at
% Austin's Institute for Geophysics (UTIG), January 2016. 
% http://www.chadagreene.com 
% 
% See also: cell2mat, nan, and cat. 

%% Input checks:

narginchk(1,1) 
assert(iscell(A),'Input error: Input must be a cell array.')

%% Perform mathematics and whatnot: 

% Append a NaN to each array inside A: 
Anan = cellfun(@(x) [x;NaN],A,'un',0);

% Columnate: 
B = cell2mat(Anan(:));

end

%                                       I                                         
%                                   +Z7$I77O8M$$$$I                               
%                                 IM$MMMZMMMMMMMMMMIO?                            
%                                ?ZMMMMMMMMMMMMMMMMMMM77777M                      
%                              ,7$MMMMMMMMMMMMMMMMMMMMMMMMMM7,,                   
%   :                          I7$MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMD+             
%   8                         7IIMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMDMM=         
% ~,M,                        :7ZMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM         
%  ~M$                         MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM:         
%  $ZNZ7?                     MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMI         
%    DOII7I                  =MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM         
%    =MO7$Z8$Z               ~MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMN        
%     DZZMMMMMM8$          :77$IOMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM        
%         OMMMMMMMD :I777IOMM87IMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMZMMM        
%       7DNMMM7MMMMD77777777MM7ZMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMN$Z7I77I7        
%         7$MO7NMMMM77777777O77$MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM=         
%          ,Z $DMMMM87I7I777777MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMN,       
%            MZNMMMM$$777777DI78MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM      
%             8MMMMMMM$7Z7$MMMMMMMMMM:The Bedmap2 Toolbox:MMMMMMMMMMMMMMMMMMM$    
%             OMMMMMMMN$MMMMMMMMMMMMMMMMM:for Matlab:MMMMMMMMMMMMMMMMMMMMMMMMM7   
%            $MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMZ   
%            DMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM8   
%            7MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMN   
%           N7MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMO+: 
%           M7ZZ~7MMMMMMMMMMMMMMMMMMMM:by Chad A. Greene:MMMMMMMMMMMMMMMMMMMMMOI= 
%             +  7NMMMMMMMMMMMMMMMMMMZMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMZ7  
%                ?$MMMMMMMMMMMMMMMMMMMZNMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM   
%                7NMMMMMMMMMMMMMMMMMD$M777MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM    
%                ~ZMMMMMMMMMMMMMMMM$77I77777ZMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM     
%                 O8MMMMMMMMMMMMM77Z$I7777777OMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM7     
%                  DIMMMMMMMMMMMD77I777777777INMMMMMMMMMMMMMMMMMMMMMMMMMMMNMM     
%                  ,  $MMMMMMMMMMZI77777777777$MMMMMMMMMMMMMMMMMMMMMMMMMMM7       
%                      ,~NNMN7O8MD+,=77777777DMMMMMMMMMMMMMMMMMMMMMMMMMMM$        
%                           ,           ~~77ZD,MMMMMMMMMMMMMMMMMMMMMMMMMD         
%                                              OMMMMMMMMMMMMMMMMMMMMMMMM~         
%                                             ~?MMMMMMMMMMMMMMMMMMMMMMOD          
%                                              NMMMMMMMMMMMMMMMMMMMMMM            
%                                            :?MMMMMMMMMMMMMMMMMMMMMM~            
%                                            ,MMMMMMMMMMMMMMMMMMMMMM+             
%                                            =MMMMNMMMMMMMMMMMMMM,                
%                                                  +MMMMO    ?           
%
% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 


% Chad A. Greene
% Jackson School of Geosciences
% University of Texas Institute for Geophysics
% The University of Texas at Austin