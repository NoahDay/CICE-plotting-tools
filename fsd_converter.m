function [processed_data] = fsd_converter(filename,input,output,raw_data)
%FSD_CONVERTER will convert AFSDN, AFSD, FSDRAD as specified by the
%input/output arguments
cd /Users/a1724548/GitHub/CICE-plotting-tools % Move up one directory
% Get the dimensions
if input == "afsdn"
    dim = 4;
elseif input == "afsd"
    dim = 3;
elseif input == "fsdrad"
    dim = 2;
else
    error("Input argument not specified.")
end
sector = "world";
NCAT = ncread(filename,"NCAT");
NFSD = ncread(filename,"NFSD");
floe_rad_c = NCAT;

Nf = numel(NFSD);
Nc = numel(NCAT);
% Set up
    % floe_rad_l = lims(1:nfsd  )
    % floe_rad_h = lims(2:nfsd+1)
    % floe_binwidth = floe_rad_h - floe_rad_l
lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
floe_binwidth = floe_rad_h - floe_rad_l;
floe_rad_c = (floe_rad_l+floe_rad_h)/2;

if ~exist('data', 'var')
    raw_data = data_format_sector(filename,input,sector,dim);
end

% Grab the concentration data
aice = data_format_sector(filename,"aice",sector,dim);
aicen = data_format_sector(filename,"aicen",sector,dim);

% Get the grid dimensions
[nx,ny] = size(aice);

% Convert the data
if input == "afsdn"
   if output == "afsdn"
       % Do nothing
        processed_data = raw_data;

   elseif output == "afsd"
       % Convert from afsdn to afsd
       for j = 1:ny
           for i = 1:nx
               work = zeros(Nf,Nc);
               for k = 1:Nf
                   for n = 1:Nc
                        work(k,n) = raw_data(i,j,k,n);
                   end
               end
               processed_data(i,j,:) = sum(work');
           end
       end

   elseif output == "fsdrad"
       % Convert from afsdn to fsdrad
       for j = 1:ny
           for i = 1:nx
               work(i,j) = 0;
               for k = 1:Nf
                   for n = 1:Nc
                       work(i,j) = work(i,j) + ...
                           raw_data(i,j,k,n)*floe_binwidth(k)*(floe_rad_c(k)/aice(i,j));
                   end
               end
           end
       end
       processed_data = work;
        
   elseif output == "fsd"
       % Convert from afsdn to areal FSD
       % First convert to afsd
       for j = 1:ny
           for i = 1:nx
               work = zeros(Nf,Nc);
               for k = 1:Nf
                   for n = 1:Nc
                        work(k,n) = raw_data(i,j,k,n);
                   end
               end
               processed_data(i,j,:) = sum(work').*floe_binwidth;
           end
       end 

   end

elseif input == "afsd"
    if output == "afsdn"
       % Not possible
        error("Cannot convert from afsd to afsdn.")

   elseif output == "afsd"
       % Do nothing
        processed_data = raw_data;

   elseif output == "fsdrad"
        % Convert from afsd to fsdrad
        for j = 1:ny
            for i = 1:nx
                work(i,j) = 0;
                for k = 1:Nf
                     work(i,j) = work(i,j) + ...
                           raw_data(i,j,k)*floe_binwidth(k)*(floe_rad_c(k)/aice(i,j));
                end
            end
        end
        processed_data = work;
                
   elseif output == "fsd"
       % Convert from afsd to areal FSD
       processed_data = input.*floe_binwidth;

    end
elseif input == "fsdrad"
    if output == "fsd"
        % Convert from fsdrad to areal FSD
         error("Cannot convert from fsdrad to fsd.")
    end
else
    error("Input argument not specified.")
end


end