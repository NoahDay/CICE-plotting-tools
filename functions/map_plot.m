function [w, a, output_data] = map_plot(data,variable,sector,grid,clims,plotting)
%Plot a worldmap
cd '/Users/noahday/GitHub/CICE-plotting-tools'
addpath functions
if ~exist('sector', 'var')
    % Set sector to world by default
    sector = "world"; 
end

if ~exist('grid', 'var')
    % Set colormap to turbo by default
    grid = "gx1"; 
end


[lat,lon,row] = grid_read(grid);
if isstring(sector)
    if sector == "world"     
        w = worldmap('world');
            axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
            setm(w, 'Origin', [-90 0 0]);
            setm(w, 'maplatlimit', [-90,-50]);
            setm(w, 'maplonlimit', [-180,180]);
            setm(w, 'meridianlabel', 'on')
            setm(w, 'parallellabel', 'off')
            setm(w, 'mlabellocation', 60);
            setm(w, 'plabellocation', 10);
            setm(w, 'mlabelparallel', -45);
            setm(w, 'mlinelimit', [-90 0]);
            setm(w, 'grid', 'on');
            %setm(w, 'frame', 'on');
            setm(w, 'labelrotation', 'on')
            pcolorm(lat,lon,data)
            land = shaperead('landareas', 'UseGeoCoords', true);
            geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
            a = colorbar;
            if ~exist('clims', 'var')
                % Set sector to world by default
                caxis(colorlims(variable));
            else
                caxis(clims)
            end
    else
        coords = sector_coords(sector);

        min_lat = min(coords(:,1));
        max_lat = max(coords(:,1));
        min_lon = min(coords(:,2));
        max_lon = max(coords(:,2));
        w = worldmap('world');
            axesm miller; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
            setm(w, 'Origin', [0 0 0]);
            setm(w, 'maplatlimit', [min_lat,max_lat]);
            setm(w, 'maplonlimit', [min_lon,max_lon]);
            setm(w, 'meridianlabel', 'on')
            setm(w, 'parallellabel', 'on')
            setm(w, 'mlabellocation', 10);
            setm(w, 'plabellocation', 10);
            setm(w, 'mlabelparallel', 0);
            setm(w, 'mlabelParallel', 'south');
            setm(w, 'grid', 'on');
            setm(w, 'frame', 'on');
            setm(w, 'labelrotation', 'on')
            pcolorm(lat,lon,data)
            land = shaperead('landareas', 'UseGeoCoords', true);
            geoshow(w, land, 'FaceColor', [0.5 0.7 0.5]);
            a = colorbar;
            a.TickLabelInterpreter = 'latex';
            a.Label.String = colorlabel(variable);
            %colormap turbo

            if ~exist('clims', 'var')
                % Set sector to world by default
                caxis(colorlims(variable));       
            else
                caxis(clims)
            end
    end
    output_data = data;
else % Transect
     [lat_out,lon_out] = lat_lon_finder(sector(1),sector(2),lat,lon);
     [len, wid] = size(lat);
     transect_data = zeros(len,wid);
     transect_data(isnan(data)) = data(isnan(data));
     transect_data(lon_out,:) = data(lon_out,:);
     if plotting == 1
         w = worldmap('world');
                axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
                setm(w, 'Origin', [-90 0 0]);
                setm(w, 'maplatlimit', [-90,-50]);
                setm(w, 'maplonlimit', [-180,180]);
                setm(w, 'meridianlabel', 'on')
                setm(w, 'parallellabel', 'off')
                setm(w, 'mlabellocation', 60);
                setm(w, 'plabellocation', 10);
                setm(w, 'mlabelparallel', -45);
                setm(w, 'mlinelimit', [-90 0]);
                setm(w, 'grid', 'on');
                %setm(w, 'frame', 'on');
                setm(w, 'labelrotation', 'on')
                pcolorm(lat,lon,transect_data)
                land = shaperead('landareas', 'UseGeoCoords', true);
                geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
                a = colorbar;
                if ~exist('clims', 'var')
                    % Set sector to world by default
                    caxis(colorlims(variable));
                else
                    caxis(clims)
                end
     else
         w = 0;
         a = 0;
     end
            output_data = transect_data;
end


end

