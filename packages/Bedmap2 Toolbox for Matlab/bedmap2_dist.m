function d = bedmap2_dist(lati_or_xi,loni_or_yi,reference)
% bedmap2_dist calculates distance in kilometers from any location(s)
% to any Bedmap2 mask type.  It may be useful for finding distances to 
% the nearest grounding line or coast line.  
% 
% This function requires Matlab's Image Processing Toolbox. 
% 
%% Syntax
% 
%  d = bedmap2_dist(lati,loni,reference) 
%  d = bedmap2_dist(xi,yi,reference) 
% 
%% Description
% 
% d = bedmap2_dist(lati,loni,reference) returns distance(s) d from geo location(s)
% lati,loni to the nearest reference mask type. Distances are in kilometers. Options 
% for reference can be
% 
%   'grounded'   distance to nearest grounded ice. 
%   'ice sheet'  distance to nearest bit of ice sheet (grounded or ice shelf).   
%   'ice shelf'  distance to nearest floating ice shelf. 
%   'open ocean' distance to nearest open ocean (returns nonzero distance on an ice shelf).  
%   'tidal'      distance to nearest ocean or ice shelf (e.g., returns zero at the center of Ross Ice Shelf). 
%
% d = bedmap2_dist(xi,yi,reference) as above, but where input locations are in polar stereographic
% meters. Input coordinates are automatically determined by the islatlon function.  
% 
%% Example
% Consider (80°S,175°W), which is somewhere in the middle of the Ross Ice Shelf.  How far is it from
% grounded ice? 
% 
%     bedmap2_dist(-80,-175,'grounded')
%     ans =
%             223.88
% 
% That's about 224 km from (80°S,175°W) to the nearest grounding line. And how far from the nearest open ocean? 
% 
%     bedmap2_dist(-80,-175,'open ocean')
%     ans =
%             175.81
% 
% And of course (80°S,175°W) is in the middle of an ice shelf so all of these will be zero: 
% 
%     bedmap2_dist(-80,-175,'ice shelf');
%     bedmap2_dist(-80,-175,'tidal');
%     bedmap2_dist(-80,-175,'ice sheet');
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
% This function was written by Chad A. Greene of the University of Texas at
% Austin's Institute for Geophyiscs (UTIG) in December 2015. 
% http://www.chadagreene.com
% 
% See also bedmap2_interp and pathdist. 

%% Error checks: 

narginchk(3,3) 
assert(license('test','image_toolbox')==1,'The bedmap2_dist function requires Matlab''s Image Processing Toolbox.')
assert(exist('islatlon.m','file')==2,'Cannot find islatlon, which is required function for bedmap2_dist. Ensure you have Antarctic Mapping Tools (AMT), AMT is up to date, and Matlab can find the AMT functions. Get the latest AMT here:http://www.mathworks.com/matlabcentral/fileexchange/47638 if necessary.') 
assert(isnumeric(lati_or_xi)==1,'Input error: coordinates must be numeric.') 
assert(isnumeric(loni_or_yi)==1,'Input error: coordinates must be numeric.') 
assert(isequal(size(lati_or_xi),size(loni_or_yi))==1,'Input error: coordinates must be equal size.') 
assert(isnumeric(reference)==0,'Input error: reference mask type must be a string.') 

%% Parse inputs: 

if islatlon(lati_or_xi,loni_or_yi) 
    [xi,yi] = ll2ps(lati_or_xi,loni_or_yi); 
else
    xi = lati_or_xi; 
    yi = loni_or_yi; 
end

%% Import data and calculate gridded distances to reference mask: 

[X,Y,icemask] = bedmap2_data('icemask',xi,yi,1600,'xy'); 

switch lower(reference)
    case 'grounded'
        mask = icemask==0; 
    case 'ice sheet'
        mask = icemask~=127; 
    case 'ice shelf'
        mask = icemask==1; 
    case 'open ocean'
        mask = icemask==127; 
    case 'tidal'
        mask = icemask~=0; 
    otherwise
        error(['Unrecognized reference mask type ',reference,'.']) 
end

dst = double(bwdist(mask)); 

%% Interpolate: 

d = interp2(X,Y,dst,xi,yi); 

end

