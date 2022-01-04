Bedmap2 datasets and associated grids (5th March 2013)
''''''''''''''''''''''''''''''''''''''''''''''''''''''

Full details of bedmap2 can be found in this paper: http://www.the-cryosphere.net/7/375/2013/tc-7-375-2013.pdf

Please cite this paper when referring to bedmap2.


These products supersede the earlier release of draft bedmap2 grids, which have the word draft in the filename. 
We have made numerous minor improvements and corrections to the draft products, and these should no longer be used.


You have a choice of four data formats: tiff, ESRI geodatabase (gdb), ASCII (txt) or binary (bin).

For each format, there is a zip file containing ten datasets:

1. bedmap2_bed 				is bed height
2. bedmap2_surface 			is surface height
3. bedmap2_thickness 			is ice thickness
4. bedmap2_icemask_grounded_and_shelves is a mask file showing the grounding line and the extent of the floating ice shelves
5. bedmap2_rockmask 			is a mask file showing rock outcrops
6. bedmap2_lakemask_vostok 		is a mask file showing the extent of the lake cavity of Lake Vostok
7. bedmap2_bed_uncertainty 		is the bed uncertainty grid shown in figure 12 of the manuscript
8. bedmap2_thickness_uncertainty_5km 	is the thickness uncertainty grid shown in figure 11 of the manuscript
9. bedmap2_data_coverage 		is a binary grid showing the distribution of ice thickness data used in the grid of ice thickness
10. gl04c_geoid_to_wgs84 		gives the values (as floating point) used to convert from heights relative to WGS84 datum to heights relative to EIGEN-GL04C geoid (to convert back to WGS84, add this grid)


Grid projection and extents
'''''''''''''''''''''''''''
Each dataset is projected in Antarctic Polar Stereographic projection, latitude of true scale -71 degrees south, datum WGS84.
All heights are in metres relative to sea level as defined by the g104c geoid.

The grid dimensions are 6667 x 6667 cells and the extent is:
Top: 3333500
Left: -3333500
Right: 3333500
Bottom: -3333500

The bedmap2 grid spacing is 1000 m.

Exceptions are the 5 km thickness uncertainty grid (1361 x 1361, top: 3402500, left: -3401500, right: 3403500, , bottom: -3402500), 
and the Lake Vostok grid (281 x 112, top: -290500, left: 1189500, right: 1470500, bottom: -402500).

Uncertainty values are in metres.
 

Data formats
''''''''''''
The tiff, ASCII and gdb data format is 16 bit signed integer for the data grids, the masks are 1 bit or 8 bit integer.
The binary format is little-endian 32-bit single-precision floating point.

To read the binary format into Matlab, use e.g.:

fid=fopen('bedmap2_bed.flt','r','l');
bed=fread(fid,[6667,6667],'float32');
fclose(fid);


