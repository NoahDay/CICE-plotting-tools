close all
clear all
addpath functions

%Example 1, short distance: 
latlon1=[-43 172]; 
latlon2=[-44 171]; 
[d1km d2km] = distance(latlon1,latlon2) 

%d1km approximately equal to d2km

%--Example 2, longer distance: 
latlon1=[-43 172]; 
latlon2=[20 -108]; 
[d1km d2km]=distance(latlon1,latlon2) 
[d1km d2km]=lldistkm(latlon1,latlon2) 

%d1km is significantly different from d2km (d2km is not able to work for longer distances).