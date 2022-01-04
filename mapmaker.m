 function map_out = mapmaker(data,lat,lon,plot_title,map_type,variable,endline)
        color_map = seaicecolormap();
        w = worldmap('world');
        if map_type == "cassini"
            x_origin = -20;
            axesm cassini;
            setm(w, 'mlabelparallel', -45);
            setm(w, 'mlabellocation', 30);
        elseif map_type == "eqdcylin"
            x_origin = 0;
            axesm eqdcylin;
            setm(w, 'mlabelparallel', -85);
            setm(w, 'mlabellocation', 60);
        else
            x_origin = -90;
            axesm eqaazim; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim. cassini
            setm(w, 'mlabelparallel', -45);
            setm(w, 'mlabellocation', 30);
        end
        setm(w, 'Origin', [x_origin 0 0]);
        setm(w, 'maplatlimit', [-90,-40]);
        setm(w, 'maplonlimit', [-180,180]);
        setm(w, 'meridianlabel', 'on')
        setm(w, 'parallellabel', 'off')
        setm(w, 'plabellocation', 10);
        setm(w, 'mlabelparallel', -45);
        setm(w, 'grid', 'on');
        setm(w, 'labelrotation', 'on')
        pcolorm(lat,lon,data)
        land = shaperead('landareas', 'UseGeoCoords', true);
        geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])
        title(plot_title, 'Units', 'normalized', 'Position', [-0.0, 0.5, 0], 'Rotation',90,'Interpreter', 'none')%title(plot_title, 'Interpreter', 'none',[0.5, -0.1, 0])
        limits = colorlims(variable);
        caxis(limits)
        if endline == 1
            colorbar
        end
    end