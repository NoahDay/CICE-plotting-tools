

%% Requirements 
% 
% 1. GET ANTARCTIC MAPPING TOOLS:
% http://www.mathworks.com/matlabcentral/fileexchange/47638

% 2. GET BEDMAP2 DATA:
% Two options for getting Bedmap2 data: You can try running this script, 
% and that might should be it. But if this download script fails, you can manually
% download the data from here: https://secure.antarctica.ac.uk/data/bedmap2/bedmap2_tiff.zip, 
% be sure to unzip it and put it somewhere Matlab can find it. 

%% Citing Bedmap2
% If the Bedmap2 Toolbox iis useful for you, please cite the following: 
% 
% Fretwell, P., et al. "Bedmap2: improved ice bed, surface and thickness datasets for Antarctica." The Cryosphere 7.1 (2013).
% http://dx.doi.org/10.5194/tc-7-375-2013
% 
% Greene, C. A., Gwyther, D. E., & Blankenship, D. D. Antarctic Mapping Tools for Matlab. 
% Computers & Geosciences. 104 (2017) pp.151-157. http://dx.doi.org/10.1016/j.cageo.2016.08.003
% 

%% Download data:
% This is a large data set (222 MB). Depending on your internet connection speed
% it may take a while to download.

if exist('bedmap2_bed.tif','file')==0 % If an unzipped tiff does not exist, look for a zipped-up version.
    if exist('bedmap2_tiff.zip','file')~=0 % If a zipped file exists, unzip it. 
        unzip('bedmap2_tiff.zip'); 
    else
        disp('Cannot find the bedmap2_tiff.zip data set. Downloading now...')
        try
                unzip('https://secure.antarctica.ac.uk/data/bedmap2/bedmap2_tiff.zip');
                disp('Bedmap2 data download complete.') 
            catch err
                error('MATLAB:Unzip', ...
                    ['There''s been a problem in trying to download Bedmap2 data.\n',...
                'This likely relates to your internet firewall.\n',...
                'The simplest solution is to download the bedmap2_tiff.zip file from \n',...
                'https://secure.antarctica.ac.uk/data/bedmap2/ and unzip it manually.\n',...
                'Be sure to place the unzipped files in your current directory and then \n',...
                'rerun this installation script.']);

        end
    end
end

try % Move contents of unzipped folder to 
    movefile('bedmap2_tiff/*');
end 

disp 'Download successful. I think you''re ready to use the Bedmap2 plugin.'

