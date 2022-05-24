function processed_data = fsd_converter(filename,input,output,raw_data,sector)
%FSD_CONVERTER will convert between AFSDN, AFSD, FSDRAD, and DAFSD as 
% specified by the input/output arguments.
% Inputs: (filename, input_var, output_var, raw data, sector)
% Types: (string, string, string, matrix, string)

% Get the dimensions
if input == "afsdn"
    dim = 4;
elseif input == "afsd"
    dim = 3;
elseif input == "fsdrad"
    dim = 2;
elseif input == "dafsd_newi" || input == "dafsd_latg" || input == "dafsd_latm" || input == "dafsd_wave" || input == "dafsd_weld" || input == "dafsd"
    dim = 3;
else
    error("Input argument not specified.")
end

if ~exist('sector','var')
    sector = "SH";
end

NCAT = ncread(filename,"NCAT");
NFSD = ncread(filename,"NFSD");
floe_rad_c = NCAT;

Nf = numel(NFSD);
Nc = numel(NCAT);

% Set up
lims = [6.65000000e-02,   5.31030847e+00,   1.42865861e+01,   2.90576686e+01, 5.24122136e+01,   8.78691405e+01,   1.39518470e+02,   2.11635752e+02, 3.08037274e+02,   4.31203059e+02,   5.81277225e+02,   7.55141047e+02, 9.45812834e+02,   1.34354446e+03,   1.82265364e+03,   2.47261361e+03,  3.35434988e+03];
floe_rad_l = [lims(1:Nf)]; % Floe radius lower bound
floe_rad_h = lims(2:Nf+1); % Floe radius higher bound
floe_binwidth = floe_rad_h - floe_rad_l;
floe_rad_c = (floe_rad_l+floe_rad_h)/2;

% Grab the concentration data
if output == "aice"
    nx = 321;
    ny = 384;
else
    aice = data_format_sector(filename,"aice",sector);
    % Get the grid dimensions
    [nx,ny] = size(aice);
end
%aicen = data_format_sector(filename,"aicen",sector,dim);



if ~exist('raw_data', 'var')
    raw_data = data_format_sector(filename,input,sector);
end



% Convert the data
if input == "afsdn"
   if output == "afsdn"
       % Do nothing
        processed_data = raw_data;

   elseif output == "afsd"
       % Convert from afsdn to afsd
       % Integrate the joint floe size and thickness distribution to obtain
       % the FSD, F(r)
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
       % Convert from f(r,h) to the probability mass function L(r,h)dr
       % where L(r,h)dr is the fraction of ice with lateral floe size
       % between r and r+dr and satisfies sum(L(r,h)) = 1
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


   elseif output == "n_fstd"
      % Calculate the number fsd, units (km^{-2})
       alpha = 0.66; % Rothrock and Thorndike (1984)
       % Get tracer array
       aicen_data = data_format_sector(filename,"aicen",sector);
        for i = 1:nx
            for j = 1:ny
                for k = 1:nf
                    for n = 1:nc
                        trcrn(i,j,k,n) = raw_data(i,j,k,n)*floe_binwidth(k)/aicen_data(i,j,n); % something x Length: (something m)
                    end
                end
            end
        end
        % Calculate nfstd
        alpha = 0.66; % Dimensionless
        floe_area_c = 4*alpha*floe_rad_c.^2; % Area: (m^2)
        for i = 1:nx
            for j = 1:ny
                for k = 1:nf
                    for n = 1:nc
                        processed_data(i,j,k,n) = trcrn(i,j,k,n)/floe_area_c(k); % Dimensionless: something / Length
                    end
                end
            end
        end

   elseif output == "n_fsd"
        alpha = 0.66; % Rothrock and Thorndike (1984)
       % Get tracer array
       aicen_data(:,:,:) = data_format_sector(filename,"aicen",sector);
       afsdn = data_format_sector(filename,"afsdn",sector);
%       error(sprintf('%d\n',size(afsdn)))
        for i = 1:nx
            for j = 1:ny
                for k = 1:Nf
                    for n = 1:Nc
                        trcrn(i,j,k,n) = raw_data(i,j,k,n).*floe_binwidth(k)*(10^(-3))./aicen_data(i,j,n); % something x Length: (something m)
                    end
                end
            end
        end
        % Calculate nfstd
        alpha = 0.66; % Dimensionless
        floe_area_c = 4*alpha*(floe_rad_c*(10^(-3))).^2; % Area: (m^2)
        for i = 1:nx
            for j = 1:ny
                for k = 1:Nf
                    for n = 1:Nc
                        nfstd(i,j,k,n) = trcrn(i,j,k,n)/floe_area_c(k); % Dimensionless: something / Length
                    end
                end
            end
        end
        % Intergrate nfstd to get nfsd

        for i = 1:nx
            for j = 1:ny
                for k = 1:Nf
                    processed_data(i,j,k) = sum(nfstd(i,j,k,:)); % number of floes per m^2
                end
            end
        end

   elseif output == "ITD"
        % Integrate w.r.t. floe size to obtain the ITD g(h)
   elseif output == "aice"
       % From AFSDN get AFSD and then aice
       for j = 1:ny
           for i = 1:nx
               work(i,j) = 0;
               for k = 1:Nf
                   for n = 1:Nc
                       work(i,j) = work(i,j) + ...
                           raw_data(i,j,k,n)*floe_binwidth(k);
                   end
               end
           end
       end
       processed_data = work;
   else
       error("No output specified")
   end

   
   % AFSD input
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
       % Convert from F(r) to the probability mass function L(r,h)dr
       % where L(r,h)dr is the fraction of ice with lateral floe size
       % between r and r+dr and satisfies sum(L(r,h)) = 1
       processed_data = input.*floe_binwidth;

    end
    
    
    % FSDRAD input
elseif input == "fsdrad"
    if output == "fsd"
        % Convert from fsdrad to areal FSD
         error("Cannot convert from fsdrad to fsd.")
    end
    
    
elseif input == "dafsd_newi" || input == "dafsd_latg" || input == "dafsd_latm" || input == "dafsd_wave" || input == "dafsd_weld" || input == "dafsd" % dafsd input
    if output == "fsdrad"
    % Apply same as afsd -> fsdrad
        for j = 1:ny
            for i = 1:nx
                work(i,j) = 0;
                for k = 1:Nf
                     work(i,j) = work(i,j) + ...
                           raw_data(i,j,k)*floe_rad_c(k)./aice(i,j);
                end
            end
        end
        processed_data = work; % normalise by ice concentration
    elseif output == "m_per_day"
        for j = 1:ny
            for i = 1:nx
                work(i,j) = 0;
                for k = 1:Nf
                     work(i,j) = work(i,j) + ...
                           raw_data(i,j,k)*floe_rad_c(k);%*floe_binwidth(k)*(floe_rad_c(k)/aice(i,j));
                end
            end
        end
        processed_data = work;
    else
        processed_data = raw_data;
    end
else
    error("Input argument not specified.")
    
end


end