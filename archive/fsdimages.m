% Plotting the FSD hsitograms along a specified transect of the CICE output
% Comparing CICE results with and without waves
clear all
close all
cases = ["8month", "1yearnowaves"]; %8month

variable = ["dafsd_wave", "dafsd_weld", "dafsd_latm", "dafsd_latg", "dafsd_newi"];
latit = 1;
% Read the header
nvars = numel(variable);
ncases = numel(cases);

%% Figure 1
figure(1)
%t1 = tiledlayout(nfsd/2);
for k = 1:ncases
    filename = strcat('cases/',cases(k),'/history/iceh.2005-08-31.nc');
    %ncdisp(filename)
    % grid
    lat = ncread('grid/global_gx1.bathy.nc','TLAT');
    lon = ncread('grid/global_gx1.bathy.nc','TLON');

    dim = 2;
    lat = rearrange_matrix(lat,37,dim);
    lon = rearrange_matrix(lon,37,dim);

    lon = [zeros(1,384);lon];
    lat = [lat(1,:); lat];
    
    fsd_cat = ncread(filename, 'NFSD');
    fsd_range = [0; fsd_cat];
    nfsd = length(fsd_cat);

   
    for j = 1:nvars
        data_read = ncread(filename, variable(j)); 
        for i = 1:nfsd % number of FSD categories
            data_2 = data_read(:,:,i);
            data = rearrange_matrix(data_2,37,dim);
            data = [data; data(end,:)];
            color_map = seaicecolormap();
            %if max(max(abs(data)))>eps
                % Mapping
                %nexttile%(i)
                latitude = [-90,-30];
                longitude = [-180,180];
                w = worldmap('world');
                axesm eqdcylin; %, eqaazim eqdazim vperspec, eqdazim flips the x-axis, and y-axis to eqaazim.
                setm(w, 'Origin', [0 0 0]);
                setm(w, 'maplatlimit', [-90,-40]);
                setm(w, 'maplonlimit', [-180,180]);
                setm(w, 'meridianlabel', 'on')
                setm(w, 'parallellabel', 'off')
                setm(w, 'mlabellocation', 60);
                setm(w, 'plabellocation', 10);
                setm(w, 'mlabelparallel', -85);
                setm(w, 'grid', 'on');
                %setm(w, 'frame', 'on');
                setm(w, 'labelrotation', 'on')
                pcolorm(lat,lon,data)
                land = shaperead('landareas', 'UseGeoCoords', true);
                geoshow(w, land, 'FaceColor', [0.5 0.7 0.5])

                single_title = sprintf('[%0.1f, %0.1f] m', fsd_range(i),fsd_range(i+1));
                double_title = sprintf('FSD category \n [%0.1f, %0.1f] m', fsd_range(i),fsd_range(i+1));
                if i == 1
                    % Add layout title
                    if variable(j) == "dafsd_wave"
                     plot_var = " waves";
                    elseif variable(j) == "dafsd_weld"
                     plot_var = " welding";
                    elseif variable(j) == "dafsd_latm"
                     plot_var = " lateral melt";
                    elseif variable(j) == "dafsd_latg"
                     plot_var = " lateral growth";
                    elseif variable(j) == "dafsd_newi"
                     plot_var = " new ice";
                    end
                    plot_title = strcat('Change in the FSD due to', plot_var, ', 31-08-2005');
                    title(plot_title,'FontSize',30)
                end
                set(gca,'FontSize',10)
                set(gcf, 'Position',  [100, 100, 1200, 200])
                text(-70, -60, '(60 N, 52 W)');
                a=colorbar;
                label_c = ylabel(a,double_title,'FontSize',12,'Rotation',90);
                label_c.Position(1) = -57;
                label_h.Position(2) = 1; % change vertical position of ylabel
                %caxis([-0.01,0.01]);

                % Save the images
                image_name = strcat(cases(k),'_',variable(j),'%d.png');
                figname = sprintf(image_name, i); 
                %text = sprintf('/Users/%s/MATLAB-Drive/MATLAB/PhD Project/CICE Plotting/frames', user);
                fname = '/Volumes/Noah SSD/MATLAB/PhD Project/CICE Plotting/frames/fsd'; 
                saveas(gcf,fullfile(fname, figname));
           % end
        end
    end
end