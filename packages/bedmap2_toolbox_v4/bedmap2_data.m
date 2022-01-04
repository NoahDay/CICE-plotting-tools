function varargout = bedmap2_data(zvar,varargin)
% bedmap2_data returns Bedmap2 data and corresponding geographic or polar
% stereographic coordinates. 
% 
% This function requires the Antarctic Mapping Tools> package found here: 
% http://www.mathworks.com/matlabcentral/fileexchange/47638
% 
%% Syntax 
% 
%  Z = bedmap2_data(var)
%  Z = bedmap2_data(...,lati,loni)
%  Z = bedmap2_data(...,lati,loni,extrakm)
%  Z = bedmap2_data(...,xi,yi,extrakm)
%  Z = bedmap2_data(...,xi,yi)
%  Z = bedmap2_data(...,'res',resolution)
%  [lat,lon] = bedmap2_data('latlon')
%  [x,y] = bedmap2_data('xy')
%  [lat,lon,Z] = bedmap2_data(...)
%  [x,y,Z] = bedmap2_data(...,'xy') 
% 
%% Description 
% 
% Z = bedmap2_data(var) returns Bedmap2 data of type var at full 1 km resolution. var options are: 
%
% * 'surface'         surface elevation of the ice (m) relative to gl04c geoid. 
% * 'surfacew'        surface elevation of the ice (m) relative to WGS84.  
% * 'bed'             bed elevation (m) relative to gl04c geoid. 
% * 'bedw'            bed elevation (m) relative to WGS84 ellipsoid. 
% * 'thickness'       ice thickness (m) 
% * 'beduncertainty'  uncertainty (m)
% * 'gl04c'           gl04c geoid (z_wgs84 = z_bedmap2 + gl04c)   
% * 'coverage'        
% * 'icemask'         0 = grounded, 1 = ice shelf, 127 = ocean
% * 'rockmask'        1s indicate rocks
% * 'vostok'          1s indicate Lake Vostok
% 
% Z = bedmap2_data(...,lati,loni) returns only enough Bedmap2 data to fully
% encompass a set of points given by geo coordinates lati,loni. 
%
% Z = bedmap2_data(...,lati,loni,extrakm) as above, but encompasses
% points lati,loni and adds a buffer of specified width extrakm in kilometers
% around all four sides of data points. 
% 
% Z = bedmap2_data(...,xi,yi) returns only enough Bedmap2 data to fully
% encompass a set of points given by polar stereographic (-71) coordinates xi,yi. 
% 
% Z = bedmap2_data(...,xi,yi,extrakm) as above, but encompasses
% points xi,yi and adds a buffer of specified width extrakm in kilometers
% around all four sides of data points. 
% 
% Z = bedmap2_data(...,'res',resolution) specifies a resolution in
% kilometers. By default, this function returns the full 1 km Bedmap2 data set. 
% 
% [lat,lon] = bedmap2_data('latlon') returns a lat,lon grid of the
% Bedmap2 data set. 
% 
% [x,y] = bedmap2_data('xy') returns a polar stereographic (-71) grid of
% the Bedmap2 data set. 
% 
% [lat,lon,Z] = bedmap2_data(...) returns geo coordinates and elevation
% or mask data Z. 
% 
% [x,y,Z] = bedmap2_data(...,'xy') returns polar stereographic
% coordinates when the 'xy' tag is included. 
% 
% 
%% Examples 
% 
% % * * * THE WHOLE BEDMAP2 DATA SET: * * * 
% 
% % Load just the surface:
% sfz = bedmap2_data('surface'); 
% 
% % Load just a lat/lon grid: 
% [lat,lon] = bedmap2_data('latlon'); 
% 
% % Load surface and lat/lon grid together: 
% [lat,lon,sfz] = bedmap2_data('surface'); 
% 
% % Load polar stereographic (-71) x,y grid: 
% [x,y] = bedmap2_data('xy'); 
% 
% % Load polar stereographic x,y grid and surface: 
% [x,y,sfz] = bedmap2_data('surface','xy'); 
% 
% 
% % * * * SUBSET BY REGION: * * * 
% 
% % Given these data points: 
% lat_gps = -71+4*rand(25,1); 
% lon_gps = -70+7*rand(25,1);
% 
% % Load only enough bed elevation data to fully contain your GPS points: 
% [lat,lon,bed] = bedmap2_data('bed',lat_gps,lon_gps); 
% 
% % Add a 100 km buffer around the GPS data points: 
% [lat,lon,bed] = bedmap2_data('bed',lat_gps,lon_gps,100); 
% 
% 
% % * * * REDUCED RESOLUTION: * * * 
% 
% % Load full data set at 5 km resolution: 
% [lat,lon,bed] = bedmap2_data('bed','res',5); 
% 
% % Regional, 2 km resolution, 75 km buffer around gps points: 
% [x,y,bed] = bedmap2_data('bed',lat_gps,lon_gps,75,'res',2,'xy'); 
% 
%% Author Info
% 
% This function and supporting documentation were written by Chad A. Greene
% of the University of Texas at Austin's Institute for Geophysics (UTIG). 
% Feel free to contact Chad by any of the means listed on his internet 
% website, http://www.chadagreene.com. 
% 
%% Change log 
% 
% March 2013: 
% * Written and uploaded to the Mathworks File Exchange site. 
% 
% December 2014: 
% * Now if user requests only one output argument, this
% function assumes the user wants the variable data.  That's because 
% [lat,lon,vardata] as an output makes sense, and vardata alone as an
% output makes sense, but lat alone does not make sense. This change was
% made for usability and improved speed. 
%
% February 2015: 
% * Switched to .tif-based data, thus eliminating the need to run
% bedmap2_install. The change to tif also brings the ability to load a
% sub-region of data, making loading and later processing data much faster. 
% 
% * Introduced regional data loading. 
% 
% * Dropped support for velocity and accumulation data. For velocity data,
% use the measures function available on the Mathworks site. If you want 
% accumulation data, let me know and maybe I'll write a function for it. Or 
% feel free to write your own accumulation function and even bill it as a 
% plugin for AMT if you'd like. 
% 
% * Dropped grounding line and coast line option, because the grounding line
% is not really a Bedmap2 product. For more accurate grounding lines, consider
% using one of the following data sets which are available as Antarctic
% Mapping Tools plugins: 
%     * measures (dInSAR-derived)
%     * icesat (laser-derived)
%     * modismoa (modis surface-slope-derived)
%     * asaid (landsat surface-slope-derived)
% If you still want grounding lines or coast lines inferred from Bedmap2
% data, type load bedmap2gl  or load bedmap2coast.
% 
% January 2016: 
% * Fixed Vostok mask.  Previously, it was actually calling rock mask 
% when user requested Lake Vostok mask data. 
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
% See also bedmap2, ll2ps, and bedmap2_interp.  

%% Initial error checks

narginchk(1,7)
assert(exist('bedmap2_bed.tif','file')==2,'Cannot find Bedmap2 data. It is possible that it''s in some folder Matlab cannot find. Another possibility is you have not updated to the .tif-based version 4 of the Bedmap2 toolbox. I do sincerely apologize for changing the underlying structure of this toolbox a couple of times now; I know that can be a real headache. But the .tif version is much faster, takes up less space on your computer, and is less of a numerical black box.')
assert(isnumeric(zvar)~=1,'The third input when using bedmap_interp must be a string like ''bed'', but it looks like you''ve tried to enter numeric values.')
assert(any(strcmpi(zvar,'gl')|strcmpi(zvar,'coast')|strcmpi(zvar,'cl')|strncmpi(zvar,'grounding',8))==0,...
    'It looks like you are trying to load a Bedmap2 coast line or grounding line. Sure, we used to do that kind of stuff around here, but not anymore. Ya see, Bedmap2 has never actually published a grounding line or a coast line. The best we can do is approximate it by converting their raster dataset into an outline. That might be fine for some applications, but if you need grounding line data, it is better to use asaid, measures, icesat, or modismoa data, which are all specifically focused on the grounding line. But if you really want Bedmap2 grounding lines or coast lines, type load bedmap2gl or load bedmap2coast.')

%% Set defaults and parse inputs

subset = false; % use whole data set (not a regional subset) by default 
extrakm = 0;    % zero buffer by default
xyout = false;  % give lat lon grid by default
res = 1000;     % full 1 km resolution by default
loadData = true;% load bedmap2 data unless xy or latlon are selected. 

if nargin>1
    if isnumeric(varargin{1}) % A clause in case varargin{1} is the 'xy' or 'res' option.
        subset = true; 
        lati_or_xi = varargin{1}; 
        loni_or_yi = varargin{2}; 

        % Are inputs georeferenced coordinates or polar stereographic?
        if islatlon(lati_or_xi,loni_or_yi)
            lati_or_xi = -abs(lati_or_xi); % Ensures the southern hemisphere.
            [xi,yi] = ll2ps(lati_or_xi,loni_or_yi); % The ll2ps function is in the Antarctic Mapping Tools package. 
        else
            xi = lati_or_xi;
            yi = loni_or_yi;  
        end

        if nargin>3
            if isscalar(varargin{3})
                extrakm = varargin{3}; % buffer to add around data, in kilometers.   
            end
        end
    end
    
    % Is the user requesting x and y outputs instead of default lat,lon grid? 
    if any(strcmpi(varargin,'xy')) 
        xyout = true; 
    end
    
    % Parse desired resolution however user may declare it: 
    tmp = strncmpi(varargin,'res',3); 
    if any(tmp)
        res = varargin{find(tmp)+1}; 
        switch res
            case {'high','hi','h','full','raw','hires','high res','1 km','1km','1',1,1000}
                res = 1000; 
                
            case {'med','medium','m','medres','midres','2 km','2km','2',2,2000}
                res = 2000; 
                
            case {'l','low','lo','low res','lores','lowres','3 km','3km','3',3,3000}
                res = 3000; 
                
            case {'4 km','4km','4',4,4000}
                res = 4000; 
                
            case {'extralow','xtralow','xtralo','lousy','xl','xlow','xlo',...
                    'xlow res','xlores','xlowres','extra low','5 km','5km','5',5,5000}
                res = 5000; 
                
            case {'10 km','10km','10',10,10000}
                res = 10e3; 
                
            otherwise
                error('Unrecognized resolution argument in bedmap2_data.')

        end
    end
end

%% Load data 

% Full Bedmap2 polar stereographic (71) grid: 
x = -3333e3:1e3:3333e3; 
y = 3333e3:-1e3:-3333e3; 

switch lower(zvar)
    case {'xy','x/y'} 
        assert(nargout==2,'If all you''re requesting is an ''xy'' grid, you must request exactly two outputs (e.g., [x,y] = bedmap2_data(''xy'') )')
        loadData = false; 
        xyout = true; 

    case {'ll','latlon','lat lon','lat/lon'} 
        assert(nargout==2,'If all you''re requesting is a ''latlon'' grid, you must request exactly two outputs (e.g., [lat,lon] = bedmap2_data(''latlon'') )')
        loadData = false; 
        
    case {'bed','bottom','bedrock','depth','bedw','bottomw','bedrockw','depthw','bedwgs','bedwigs','bedwgs84',...
            'bed wgs','bed wgs 84'}
        filename = 'bedmap2_bed.tif'; 
        
    case {'surface','sfc','surf','sfz','surfacew','sfcw','surfw','sfzw','surfacewgs','sfcwgs','surfwgs','sfzwgs',...
            'surfacewgs84','sfcwgs84','surfwgs84','sfzwgs84'}
        filename = 'bedmap2_surface.tif'; 

    case {'thickness','thick','th','thck'}
        filename = 'bedmap2_thickness.tif'; 
            
    case {'coverage','cov'}
        filename = 'bedmap2_coverage.tif';         
            
    case {'icemask','ice'} 
        filename = 'bedmap2_icemask_grounded_and_shelves.tif'; 
            
    case {'rockmask','rock'} 
        filename = 'bedmap2_rockmask.tif'; 
         
    case {'beduncertainty','bed uncertainty','bedunc','beduncadunk'} 
        filename = 'bedmap2_grounded_bed_uncertainty.tif'; 
            
    case {'lakemask','lake','vostok','vostokmask'} 
        x = 1190000:1000:1470000;
        y = -291000:-1000:-402000; 
        filename = 'bedmap2_lakemask_vostok.tif'; 
        
    case {'geoid','gl04c','wgs','wgs84'}
        filename = 'gl04c_geiod_to_WGS84.tif'; 
        
    otherwise 
        error(['Unrecognized dataset ""',zvar,'".']) 
end

if subset
    % A crude manual fix for when a single xi,yi lies between pixels: 
    if isscalar(xi)
          extrakm = max([extrakm 1]); 
    end
    
    % Get xlimits (xl) and ylimits (yl) of input coordinates + buffer:
    xl = [min(xi(:))-extrakm*1000 max(xi(:))+extrakm*1000];
    yl = [min(yi(:))-extrakm*1000 max(yi(:))+extrakm*1000];

    % Region of rows and columns of pixels to read: 
    ri=find((y>=yl(1))&(y<=yl(2)));
    ci=find((x>=xl(1))&(x<=xl(2)));

else
    ri = [1 6667]; 
    ci = [1 6667]; 
end

if loadData
    
    % Load data: 
    Z = double(imread(filename,'PixelRegion',{[min(ri) res/1000 max(ri)] [min(ci) res/1000 max(ci)]}));

    % Take care of NaNs in the tif (which are specified where values equal 32767): 
    Z(abs(Z)>32766)=NaN; 
    
    % Convert from gl04c geoid to WGS84 ellipsoid if requested by user: 
    switch lower(zvar)
        case {'bedw','bottomw','bedrockw','depthw','bedwgs','bedwigs','bedwgs84',...
                'bed wgs','bed wgs 84','surfacew','sfcw','surfw','sfzw','surfacewgs','sfcwgs','surfwgs','sfzwgs',...
                'surfacewgs84','sfcwgs84','surfwgs84','sfzwgs84'}

                Z = Z + double(imread('gl04c_geiod_to_WGS84.tif','PixelRegion',{[min(ri) res/1000 max(ri)] [min(ci) res/1000 max(ci)]}));
        otherwise
            % carry on
    end
end


%% Arrange outputs in order requested by user: 

switch nargout
    case {0,1}
        varargout{1} = Z; 
        
    case 2
        [X,Y] = meshgrid(x(min(ci):res/1000:max(ci)),y(min(ri):res/1000:max(ri)));
        
        if xyout
            varargout{1} = X; 
            varargout{2} = Y; 
        else
            [varargout{1},varargout{2}] = ps2ll(X,Y); 
        end
        
    case 3        
        [X,Y] = meshgrid(x(min(ci):res/1000:max(ci)),y(min(ri):res/1000:max(ri)));
        
        if xyout
            varargout{1} = X; 
            varargout{2} = Y;
        else
            [varargout{1},varargout{2}] = ps2ll(X,Y); 
        end
        varargout{3} = Z; 
        
    otherwise
        error('Too many outputs.') 

end
end



