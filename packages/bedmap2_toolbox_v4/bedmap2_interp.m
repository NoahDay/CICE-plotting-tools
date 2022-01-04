function vari = bedmap2_interp(lati_or_xi,loni_or_yi,zvar,method)
% bedmap2_interp returns Bedmap2 data, interpolated to any lat,lon point,
% line, or grid. This function is part of the Bedmap2 Toolbox for Matlab, 
% which is a plugin for the Antarctic Mapping Tools package. 
% 
%% Syntax
% 
% vari = bedmap2_interp(lati,loni,var); 
% vari = bedmap2_interp(xi,yi,var); 
% vari = bedmap2_interp(...,method); 
% 
%% Description
% 
% vari = bedmap2_interp(lati,loni,var) returns interpolated values of the
% variable var at locations given by geographic coordinates lati and loni.
% lati and loni may be scalar, vector, or matrix, but their dimensions must
% match. Input var can be any of the following: 
% 
%     'surface'         surface elevation of the ice (m) relative to gl04c geoid. 
%     'surfacew'        surface elevation of the ice (m) relative to WGS84.  
%     'bed'             bed elevation (m) relative to gl04c geoid. 
%     'bedw'            bed elevation (m) relative to WGS84 ellipsoid. 
%     'beduncertainty'  uncertainty (m)
%     'thickness'       ice thickness (m) 
%     'gl04c'           gl04c geoid (z_wgs84 = z_bedmap2 + gl04c)   
%     'coverage'        
%     'icemask'         0 = grounded, 1 = ice shelf, 127 = ocean 
%     'rockmask'        indicates rocks
%     'vostok'          indicates Vostoks 
% 
% vari = bedmap2_interp(xi,yi,var) accepts inputs as polar stereographic
% (true lat -71) x,y locations in meters. If any the first input to bedmap2_interp 
% contains *any* elements whose absolute value exceeds 90, polar stereographic units
% are assumed. 
% 
% vari = bedmap2_interp(...,method) specifies an interpolation
% method, which can be 'nearest','linear','spline', or'cubic'. Default 
% interpolation method is 'linear'. 
% 
%% Examples
%
% Example 1: Find the thickness of the ice at Byrd ice core in West Antarctica: 
%
% [byrdlat,byrdlon] = coreloc('byrd');
% ByrdCoreDepth = bedmap2_interp(byrdlat,byrdlon,'thickness')
% ByrdCoreDepth =
%        2147.78
% 
% Example 2: Regrid bed data to a grid defined by mylat and mylon: 
% 
% bed_mygrid = bedmap2_interp(mylat,mylon,'bed');
% 
% 
% Example 3: Cubic interpolation along a line, roughly 10 m spacing:
% 
% lats = -80:-1/11111.1:-81; 
% lons = -90*ones(size(lats)); 
% bedCubic = bedmap2_interp(lats,lons,'bed','cubic');
% plot(lats,bedCubic)
%
%% Author Info and revision history
% 
% This script and supporting documentation were written by Chad A. Greene 
% of the University of Texas at Austin's Institute for Geophysics (UTIG). 
% http://www.chadagreene.com.
% 
% June 2013: First version written for inclusion in the Bedmap2 Toolbox 
% for Matlab. 
% 
% September 2014: Fixed a half-pixel (500 m) offset. If you downloaded a 
% this script from the Mathworks website between September 2014 and
% February 2015, you likely had a version containing the error. If you got 
% the toolbox directly from Chad in that time, you probably had a fixed
% version. 
% 
% February 2015: Switched from a reformatted .mat version of the data to
% now loading the necessary subregion of the .tif data. This is much
% faster, especially when lati and loni are contained in a small region of
% Antarctica. Also dropped support for velocity data because the
% measures_interp function (available as a plugin for Antarctic Mapping
% Tools) is more powerful and dropped support for Arthern's accumulation 
% data because in the spirit of having specific plugins for Antarctic Mapping
% Tools, it does not make much sense to include accumulation with the
% Bedmap2 data set. If you miss accumulation data let me know and I might
% bring it back in some other form. This version also now accepts polar
% stereographic units as inputs, and no longer requires running
% bedmap2_install. 
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
% See also bedmap2, measures_interp, bedmap2_profile, and bedmap2_data. 

%% Error checks: 

narginchk(3,4)
assert(isnumeric(lati_or_xi)==1,'The first input when using bedmap2_interp must be numeric latitudes.') 
assert(isnumeric(zvar)~=1,'The third input when using bedmap_interp must be a string like ''bed'', but it looks like you''ve tried to enter numeric values.')
assert(numel(lati_or_xi)==numel(loni_or_yi),'Dimensions of lati and loni must match.') 

%% Set defaults and parse inputs: 

interpMethod = 'linear'; 
if nargin>3
    interpMethod = method; 
end

%% Convert inputs to polar stereographic units if necessary: 
% If any values in input 1 have an absolute value exceeding 90, input units are interpreted
% as polar stereographic (true lat -71) units. 

if islatlon(lati_or_xi,loni_or_yi)
    [xi,yi] = ll2ps(lati_or_xi,loni_or_yi); % The ll2ps function is in the Antarctic Mapping Tools package. 
else
    xi = lati_or_xi;
    yi = loni_or_yi;    
end

%% Load data via bedmap2_data: 

[X,Y,Z] = bedmap2_data(zvar,lati_or_xi,loni_or_yi,2,'xy'); 

%% Interpolate: 

vari = interp2(X,Y,Z,xi,yi,interpMethod);

end

function tf = islatlon(lat,lon)
% islatlon determines whether lat,lon is likely to represent geographical
% coordinates. 
% 
%% Syntax
% 
% tf = islatlon(lat,lon) returns true if all values in lat are numeric
% between -90 and 90 inclusive, and all values in lon are numeric between 
% -180 and 360 inclusive. 
% 
%% Example 1: A single location
% 
% islatlon(110,30)
%    = 0
% 
% because 110 is outside the bounds of latitude values. 
% 
%% Example 2: A grid
% 
% [lon,lat] = meshgrid(-180:180,90:-1:-90); 
% 
% islatlon(lat,lon)
%    = 1 
% 
% because all values in lat are between -90 and 90, and all values in lon
% are between -180 and 360.  What if it's really, really close? What if
% just one value violates these rules? 
% 
% lon(1) = -180.002; 
% 
% islatlon(lat,lon)
%    = 0
% 
%% Author Info
% This function was written by Chad A. Greene of the University of Texas at
% Austin's Institute for Geophysics (UTIG). http://www.chadagreene.com. 
% March 30, 2015. 
% 
% See also wrapTo180, wrapTo360, projfwd, and projinv.  

% Make sure there are two inputs: 
narginchk(2,2)

% Set default output: 
tf = true; 

%% If *any* inputs don't look like lat,lon, assume none of them are lat,lon. 

if ~isnumeric(lat)
    tf = false; 
    return
end

if ~isnumeric(lon)
    tf = false; 
    return
end
if any(abs(lat(:))>90)
    tf = false; 
    return
end

if any(lon(:)>360)
    tf = false; 
    return
end    

if any(lon(:)<-180)
    tf = false; 
end

end