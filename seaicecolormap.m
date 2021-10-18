function color_map = seaicecolormap()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Colormap
%                 % Nathan's colormap
        i = 1;
                CMap1 = [1, 1, 1]; % White
                CMap2 = [0, 0, 1]; % Blue
                CMap3 = [0.0588, 1, 1]; % Teal
                CMap4 = [1, 1, 0.0667]; %Yellow
                CMap5 = [1, 0, 0]; % Red
        
                num_extra_colours = 3;
        
                n = 256; %how many colours are in the colour map
        
                Colour_1 = CMap2(1); % Value for red
                Colour_2 = CMap2(2); % Value for green
                Colour_3 = CMap2(3); % Value for blue
        
                dColour1_1 = (CMap2(1)-CMap3(1)); % difference between the first, and last colour for red
                dColour1_2 = (CMap2(2)-CMap3(2)); % difference between the first, and last colour for green
                dColour1_3 = (CMap2(3)-CMap3(3)); % difference between the first, and last colour for blue
        
                Colour_4 = CMap3(1);
                Colour_5 = CMap3(2);
                Colour_6 = CMap3(3);
        
                dColour2_1 = (CMap3(1)-CMap4(1)); % difference between the first, and last colour for red
                dColour2_2 = (CMap3(2)-CMap4(2)); % difference between the first, and last colour for green
                dColour2_3 = (CMap3(3)-CMap4(3)); % difference between the first, and last colour for blue
        
                Colour_7 = CMap4(1);
                Colour_8 = CMap4(2);
                Colour_9 = CMap4(3);
        
                dColour3_1 = (CMap4(1)-CMap5(1)); % difference between the first, and last colour for red
                dColour3_2 = (CMap4(2)-CMap5(2)); % difference between the first, and last colour for green
                dColour3_3 = (CMap4(3)-CMap5(3)); % difference between the first, and last colour for blue
        
        
                for ii = i:n  % basic loop shit
                    if ii == 1
                        Mapinc = CMap1;
                        CMap1(ii,:) = Mapinc;
                    elseif ii >= 2 & ii <= n/num_extra_colours
                        ci = ii-1;
                        cm = (n/num_extra_colours)-1;
                        Mapinc = [Colour_1-dColour1_1*(ci/cm), Colour_2-(dColour1_2*(ci/cm)), Colour_3-(dColour1_3*(ci/cm))]; % each iteration of a colourmap
                        CMap1(ii,:) = Mapinc;
                    elseif ii > n/num_extra_colours & ii <= (n/num_extra_colours)*2
                        ci = ii-(n/num_extra_colours);
                        cm = (n/num_extra_colours)-1;
                        Mapinc = [Colour_4-dColour2_1*(ci/cm), Colour_5-(dColour2_2*(ci/cm)), Colour_6-(dColour2_3*(ci/cm))]; % each iteration of a colourmap
                        CMap1(ii,:) = Mapinc;
                    elseif ii > (n/num_extra_colours)*2
                        ci = ii-(n/num_extra_colours)*2;
                        cm = (n/num_extra_colours);
                        Mapinc = [Colour_7-dColour3_1*(ci/cm), Colour_8-(dColour3_2*(ci/cm)), Colour_9-(dColour3_3*(ci/cm))]; % each iteration of a colourmap
                        CMap1(ii,:) = Mapinc;
                    end
                end
        
                Map_Correction1 = find(CMap1 > 1);
                CMap1(Map_Correction1) = 1;
                Map_Correction2 = find(CMap1 < 0);
                CMap1(Map_Correction2) = 0;
        
        
                colormap (CMap1)  % sets a universal colormap, () are needed for a colour map that isn't pre-programmed.
color_map = "SeaIceColorMap";
end

