function [hice,hbed,hwater] = bedmap2_profile(lat,lon,varargin)
% bedmap2_profile plots a 2D profile of ice, water, and rock elevations
% along any lat,lon path in the Antarctic.
% 
%% Syntax
% 
%  bedmap2_profile
%  bedmap2_profile(lat,lon)
%  bedmap2_profile(...,'horiz',HorizontalAxisValues)
%  bedmap2_profile(...,PatchProperty,PatchValue)
%  bedmap2_profile(...,'inset',InsetLocation)
%  bedmap2_profile(...,'wgs84')
%  [hice,hbed,hwater] = bedmap2_profile(...)
% 
%% Description
% 
% bedmap2_profile called without any arguments opens a user interface on a
% current map. While the interface is open, click on the map to define a profile
% path. Undo points by hit *Backspace*. When you're satisfied with a path
% you've drawn, hit *Enter* to create a profile. To quit the user interface
% without creating a profile, hit *Esc*.  
% 
% bedmap2_profile(lat,lon) plots a side-view profile along a path given by
% lat,lon.  lat and lon must be 1D arrays of equal length. If only two
% points are entered, an equally-spaced 1000-point line is created between those points. 
%
% bedmap2_profile(...,'horiz',HorizontalAxisValues) specifies horizontal axis
% values where HorizontalAxisValues is a 1D monotonically-increasing or decreasing
% array of dimensions corresponding to lat and lon. By default,
% HorizontalAxisValues are calculated as pathdist(lat,lon,'kilometers') if 
% you have Matlab's Mapping Toolbox or the sum of linear distances between points
% if you do not have the Mapping Toolbox license. If you prefer to
% plot profiles with respect to some other values such as latitude of a north/south
% transect, use bedmap2_profile(lat,lon,'horiz',lat). 
%
% bedmap2_profile(...,PatchProperty,PatchValue) specifies edge line width, face
% color, and edge color of ice, water, or bed. The following properties may
% be specified: 
% 
%     * 'iceface',ColorSpec
%     * 'iceedge',ColorSpec
%     * 'iceedgewidth',LineWidth
%     * 'waterface',ColorSpec
%     * 'wateredge',ColorSpec
%     * 'wateredgewidth',LineWidth
%     * 'bedface',ColorSpec
%     * 'bededge',ColorSpec
%     * 'bededgewidth',LineWidth
%     * 'sky',ColorSpec
% 
% bedmap2_profile(...,'wgs84') plots profiles relative to the WGS84
% ellipsoid. (Profiles are plotted relative to the GL04c geoid by default.) 
% Note: Ice surface, ice thickness, and bed elevations are plotted correctly
% when using the 'wgs84' tag; however, ocean surfaces will appear at zero, 
% which could be innacurate by 60 meters or so.
% 
% bedmap2_profile(...,'inset',InsetLocation) places an inset map on the
% plot for context. Inset map placement may be specified as
% 
%     * 'NorthEast' or 'ne' for top right corner (default)
%     * 'NorthWest' or 'nw' for top left corner
%     * 'SouthWest' or 'sw' for lower left corner
%     * 'SouthEast' or 'se' for lower right corner
% 
% [hice,hbed,hwater] = bedmap2_profile(...) returns handles of ice, bed,
% and water patch objects. 
%
%% Requirements 
% 
% This function requires Matlab, Antarctic Mapping Tools, and the Bedmap2
% Toolbox for Matlab. Antarctic Mapping Tools and the Bedmap2 Toolbox 
% are available on the Mathworks File Exchange site. 
% 
%% Example 1: Profile down a glacier
% % Start by creating a crudely approximate flowline down Pine Island Glacier: 
% 
% lat = linspace(-75.447,-75.4366,50)';
% lon = linspace(-94.15,-97.951,50)';
% lat(end+1:end+50) = linspace(-75.4366,-74.9547,50);
% lon(end+1:end+50) = linspace(-97.951,-101.857,50);
% 
% bedmap2_profile(lat,lon)
% 
%% Example 2: Customize appearance
% Using the lat,lon data we generated in Example 1, we now create
% a profile of green ice on bedrock outlined by thick orange lines. And
% instead of plotting the profile as a function of distance along the path,
% we plot as a function of latitude. I'm using my rgb function to generate
% RGB triplets for colors by their names, but you can enter any RGB values 
% you'd like, or Matlab color strings like 'red' or 'w'. 
%  
% bedmap2_profile(lat,lon,...
%     'horiz',lat,...            % relative to latitude
%     'iceface',rgb('green'),... % green ice
%     'iceedge','none',...       % no ice edge color
%     'bededge',rgb('orange'),...% orange bed edge
%     'bededgewidth',8)          % thick bed edge color
% xlabel('latitude (degrees)')
% 
%% Author Info
% 
% This script was written by Chad A. Greene of the University of Texas
% at Austin's Institute for Geophysics (UTIG). February 2015. 
% http://www.chadagreene.com
% 
% July 2015: Fixed a bug whereby if no point along a profile contained surface 
% values (i.e., all ocean), the function the function previously returned an error. 
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
% See also bedmap2_interp, bedmap2_data, and pathdist. 

%% Input checks: 

if nargin>0
    assert(isvector(lat)==1,'Profile plotting requires that input latitudes in the form of a vector.'); 
    assert(isvector(lon)==1,'Profile plotting requires that input longitudes in the form of a vector.'); 
    assert(isscalar(lat)==0,'You must enter more than one point for a profile.') 
    assert(numel(lat)==numel(lon),'Input lat vector must be the same size as input lon.');
end

%% User interface if a map is already open and no inputs are defined: 

GeoMapIsOpen = false; 
psMapIsOpen = false; 
if nargin==0
    if license('test','map_toolbox') 
        if ismap(gca)
            GeoMapIsOpen = true; 
            lat = []; 
            lon = []; 
            hline = []; 
        end
    end

    if ~GeoMapIsOpen && sum(axis==[0 1 0 1])~=4
        psMapIsOpen = true; 
            xx = []; 
            yy = []; 
            hline = []; 
    end
    
    % If no map is open, open one: 
    if ~any([GeoMapIsOpen psMapIsOpen])
        bedmap2 patchshelves
        bedmap2 patchgl
        if license('test','map_toolbox') 
            GeoMapIsOpen = true; 
            lat = []; 
            lon = []; 
            hline = []; 
        else
            psMapIsOpen = true; 
            xx = []; 
            yy = []; 
            hline = []; 
        end
    end
    

    % Temporarily change figure title: 
    InitialNumberTitle = get(gcf,'NumberTitle'); 
    InitialName = get(gcf,'Name'); 
    set(gcf,'Name','Click on map to define profile:','NumberTitle','off'); 
    
    % Start GUI: 
    while 1==1
        w = waitforbuttonpress; 
        if w==1
            key = get(gcf, 'CurrentCharacter'); 

            switch key

                case 27 % 27 is the escape key, which indicates everything should be returned to initial state
                    set(gcf,'Name',InitialName,'NumberTitle',InitialNumberTitle); 
                    try; delete(hline);end
                    return 

                case 13 % The return key is 13, which means the user is happy
                    set(gcf,'Name',InitialName,'NumberTitle',InitialNumberTitle); 
                    figure
                    break

                case 8 % 8 is backspace
                    try
                        if GeoMapIsOpen
                            lat(end) = []; 
                            lon(end) = []; 
                        else
                            xx(end) = []; 
                            yy(end) = []; 
                        end
                            delete(hline(end))                            
                    end
                otherwise
                    % wait for a better command
            end
        else

            pt = get(gca,'CurrentPoint');
            if GeoMapIsOpen
                [lat(length(lat)+1),lon(length(lon)+1)] = minvtran(pt(1,1),pt(1,2));
            else
                xx(length(xx)+1) = pt(1,1); 
                yy(length(yy)+1) = pt(1,2); 
            end
        end

        try; delete(hline); end
        if GeoMapIsOpen
            hline = plotm(lat,lon,'g','linewidth',2); 
        else
            plot(xx,yy,'g','linewidth',2); 
        end
        drawnow
    end
    
end

%% Reformat lats and lons if there are not many of them: 

% Convert polar stereographic clicks to geo coordinates: 
if psMapIsOpen
    [lat,lon] = ps2ll(xx,yy); 
end

% If only two input points, turn them into a 1000 point line: 
if length(lat)==2; 
    [xi,yi] = ll2ps(lat,lon); 
    xi = linspace(xi(1),xi(2),1000); 
    yi = linspace(yi(1),yi(2),1000); 
    [lat,lon] = ps2ll(xi,yi); 
    GeoMapIsOpen = false; 
    psMapIsOpen = false; 
end

% If points are from mouse clicks, densify them, 200 points between each click: 
if any([GeoMapIsOpen  psMapIsOpen])
    [xi,yi] = ll2ps(lat,lon); 
    tmpx = [];
    tmpy = []; 
    for n = 1:length(xi)-1
        tmpx(length(tmpx)+1:length(tmpx)+200) = linspace(xi(n),xi(n+1),200); 
        tmpy(length(tmpy)+1:length(tmpy)+200) = linspace(yi(n),yi(n+1),200); 
    end
    [lat,lon] = ps2ll(tmpx,tmpy); 
end
    
% Columnate: 
lat = lat(:); 
lon = lon(:); 

%% Set defaults: 

% Colors: 
iceface = [.72 .79 .89];  
iceedge = [0 .01 .36];  
waterface = [0.0549 0.5294 0.8000]; % rgb('water blue')
wateredge = 'none'; 
bedface = [0.2039 0.1098 0.0078];   % rgb('dark brown') 
bededge = [0.1137 0.0078      0];   % rgb('very dark brown')
sky = 'w'; 

% Line widths: 
iceedgewidth = 1; 
wateredgewidth = 1; 
bededgewidth = 1; 

% Units of horizontal axis: 
distax = true; % use horizontal distance axis by by default.
if license('test','map_toolbox')
    horizax = pathdist(lat,lon,'kilometers'); % pathdist is in the Antarctic Mapping Tools package and requires Matlab's Mapping Toolbox. 
else
    % Using cumulative sum of linear distances between points is pretty
    % close usually.  I find about 2% disagreement with the great circle
    % route calculate by pathdist. If user does not have a Mapping Toolbox
    % license, cumsum should suffice: 
    [x,y] = ll2ps(lat,lon); 
    horizax = [0;cumsum(hypot(diff(x(:)),diff(y(:))))]/1000; 
end

% Inset: 
inset = false; 
location = 'ne'; % "northeast" is top right corner.  

% WGS84: 
wgs84 = false; 

%% Parse user inputs: 

if nargin>2
    
    % Horizontal Axis Values: 
    tmp = strncmpi(varargin,'horizontalaxis',3); 
    if any(tmp)
        distax = false; 
        horizax = varargin{find(tmp)+1}; 
        
        % Columnate horizontal axis to ensure consistent behavior: 
        horizax = horizax(:); 
        assert(isnumeric(horizax)==1,'horizontal axis values must be numeric.')
        assert(numel(horizax)==numel(lat),'If you enter values for a horizontal axis, they must correspond to input lat,lon values. It looks like you have entered a horizontal axis array that does not match the size of lat and lon.') 
    end
    
    % Inset: 
    tmp = strncmpi(varargin,'inset',3);
    if any(tmp)
        inset = true; 
        assert(license('test','map_toolbox')==1,'You must have a Mapping Toolbox license to include an inset map.')
        try
            location = varargin{find(tmp)+1}; 
        end
    end
    
    % WGS84 ellipsoid or GL04c geoid? (Bedmap2's GL04c geoid by default):
    if any(strncmpi(varargin,'wgs84',3))
        wgs84 = true; 
    end
    
    % All other inputs are colors or Line widths:
    tmp = strcmpi(varargin,'iceface');
    if any(tmp)
        iceface = varargin{find(tmp)+1}; 
    end
    
    tmp = strcmpi(varargin,'bedface');
    if any(tmp)
        bedface = varargin{find(tmp)+1}; 
    end
    
    tmp = strcmpi(varargin,'waterface');
    if any(tmp)
        waterface = varargin{find(tmp)+1}; 
    end
    
    tmp = strcmpi(varargin,'iceedge');
    if any(tmp)
        iceedge = varargin{find(tmp)+1}; 
    end
    
    tmp = strcmpi(varargin,'bededge');
    if any(tmp)
        bededge = varargin{find(tmp)+1}; 
    end
    
    tmp = strcmpi(varargin,'wateredge');
    if any(tmp)
        wateredge = varargin{find(tmp)+1}; 
    end
    
    tmp = strncmpi(varargin,'iceedgewidth',8);
    if any(tmp)
        iceedgewidth = varargin{find(tmp)+1}; 
        assert(isscalar(iceedgewidth)==1,'iceedgewidth value must be scalar.')
    end
    
    tmp = strncmpi(varargin,'bededgewidth',8);
    if any(tmp)
        bededgewidth = varargin{find(tmp)+1}; 
        assert(isscalar(bededgewidth)==1,'iceedgewidth value must be scalar.')
    end
    
    tmp = strncmpi(varargin,'wateredgewidth',10);
    if any(tmp)
        wateredgewidth = varargin{find(tmp)+1}; 
        assert(isscalar(wateredgewidth)==1,'iceedgewidth value must be scalar.')
    end
    
    tmp = strcmpi(varargin,'sky');
    if any(tmp)
        sky = varargin{find(tmp)+1}; 
    end
    
end

%% Get data

% Get Bedmap2 data: 
if wgs84
    sfz = bedmap2_interp(lat,lon,'surfacew','cubic'); 
    bed = bedmap2_interp(lat,lon,'bedw','cubic'); 
else
    sfz = bedmap2_interp(lat,lon,'surface','cubic'); 
    bed = bedmap2_interp(lat,lon,'bed','cubic'); 
end
thck = bedmap2_interp(lat,lon,'thickness','cubic'); 
sfz(isnan(sfz))=0; 
thck(isnan(thck)) = 0; 
icebottom = sfz-thck;

% Indices of non-ocean data: 
% nreal = isfinite(icebottom);

% Some vertical limits: 
maxsfz = max(sfz(isfinite(sfz))); 
if isempty(maxsfz)
    maxsfz = 0; 
end
minbed = min(bed(isfinite(bed))); 
padding = (maxsfz-minbed)/20; 

%% Generate plot 

% Draw water:
hwater=fill([horizax(1);horizax(end);horizax(end);horizax(1)],[0;0;minbed-padding;minbed-padding],waterface);
set(hwater,'edgecolor',wateredge,'linewidth',wateredgewidth)
hold on;

% Draw bed:
realbed = find(isfinite(bed)); 
hbed = fill([horizax(realbed);horizax(realbed(length(realbed)));horizax(realbed(1))],[bed(realbed);minbed-padding;minbed-padding],bedface);
set(hbed,'edgecolor',bededge,'linewidth',bededgewidth)

% Draw ice:
hice = patch([horizax;flipud(horizax)],[sfz;flipud(icebottom)],iceface); 
set(hice,'edgecolor',iceedge,'linewidth',1);

% Format axes:
axis tight; box off;
ylabel('elevation (m)')

% Only label x axis if it's a distance axis: 
if distax
    xlabel('along-track distance (km)')
end

% Sky color
set(gca,'color',sky)
set(gcf,'color',sky)

%% Create inset map: 

if inset
    load AMTdata.mat % contains crude grounding line and ice shelf outlines, approximated and downsampled significantly from Bedmap2
    
    gcah = gca; 
    gp = get(gca,'position');  
    
    insetsize = .25; % a quarter width of 
    insetwidth = insetsize*gp(3); 
    insetheight = insetsize*gp(4); 
    
    switch lower(location)
        case {'southwest','sw'}
            insetx = gp(1);
            insety = gp(2);   
            
        case {'northwest','nw'}
            insetx = gp(1); 
            insety = gp(2) + gp(4) - insetheight; 
            
        case {'southeast','se'}
            insetx = gp(1) + gp(3) - insetwidth;
            insety = gp(2);            
            
        otherwise % NorthEast by default:
            insetx = gp(1) + gp(3) - insetwidth;
            insety = gp(2) + gp(4) - insetheight; 
    end
    

    axes('position',[insetx insety insetwidth insetheight],'tag','insetmap');
  
        antmap('oceancolor','white')
        for k = 1:30
            patchm(slat{k},slon{k},[.8 .8 .8])
        end

        for k = 1:30
            patchm(glat{k},glon{k},[.5 .5 .5])

        end
    
    % Plot track:     
    plotm(lat,lon,'r-','linewidth',1);
    
    axes(gcah); % returns to original axes
    uistack(gcah,'down');
end


%% Clean Up

if nargout==0
    clear hice
end
