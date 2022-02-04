function [h,lat,lon] = outlineashelf(shelf,varargin)
% outlineashelf outlines the ice shelf 'shelf' on a map of Antarctica.
% 
% Now requires the Antarctic Mapping Tools package found here: 
% http://www.mathworks.com/matlabcentral/fileexchange/47638
%% Syntax 
% 
%  outlineashelf('ShelfName') 
%  outlineashelf(...,'LineProperty',LinePropertyValue)
%  h = outlineashelf(...)
%  [h,lat,lon] = outlineashelf(...)
% 
%% Description
% 
% outlineashelf('ShelfName') outlines the shelf declared by 'ShelfName'. Ice 
% shelf choices are limited to the choices below. Some ice shelves are plotted 
% together with adjacent ice shelves. 'ShelfName' can be: 
% 
%     'abbot' Abbot Ice Shelf
%     'amery' Amery Ice Shelf
%     'bach' Bach Ice Shelf
%     'baudouin' Baudouin, Borchgrevink, & Lazarev Ice Shelves
%     'conger' Conger Ice Shelf
%     'cook' Cook Ice Shelf
%     'cosgrove' Cosgrove Ice Shelf
%     'crosson' Crosson & Dotson Ice Shelf
%     'dibble' Ice Shelf
%     'dotson' Crosson & Dotson Ice Shelf
%     'drygalski' Drygalski Ice Tongue
%     'fimbul' Fimbul (Fimbulisen), Atka, & Ekstrom Ice Shelf
%     'fris' Filchner-Ronne Ice Shelf
%     'george vi' George VI Ice Shelf
%     'getz' Getz Ice Shelf
%     'holmes' Holmes Glacier
%     'land' Land & Nickerson Ice Shelf
%     'larsen b' Larsen B Ice Shelf
%     'larsen c' Larsen C & Larsen D Ice Shelf
%     'larsen d' Larsen C & Larsen D Ice Shelf
%     'larsen e' Larsen E Ice Shelf
%     'mariner' Mariner Glacier Tongue
%     'mertz' Mertz Glacier Tongue
%     'muis' Moscow University Ice Shelf
%     'nansen' Nansen Ice Shelf
%     'nickerson' Nickerson & Land Ice Shelf
%     'ninnis' Ninnis Glacier Tongue
%     'nivl' Nivl (Nivilisen) Ice Shelf
%     'pig' Pine Island Glacier shelf
%     'prince harald' Prince Harald Ice Shelf
%     'publications' Publications Ice Shelf
%     'quar' Quarisen
%     'rayner' Rayner & Thyer Ice Shelf
%     'rennick' Rennick Bay Ice Shelf
%     'riiser' Riiser-larsen, Stancomb, & Brunt Ice Shelf
%     'ross' Ross Ice Shelf
%     'shackleton' Shackleton, Tracy/Tremenchus Ice Shelf
%     'stange' Stange Ice Shelf
%     'sulzberger' Sulzberger Ice Shelf
%     'swinburne' Swinburne Ice Shelf
%     'thwaites' Thwaites Glacier Tongue
%     'totten' Totten Glacier
%     'tracy' Tracy, Tremenchus, & Shackleton Ice Shelf
%     'venable' Venable Ice Shelf
%     'vigrid' Vigridisen or Vigrid Ice Shelf
%     'vincennes' Vincennes Bay Ice Shelf
%     'west' West Ice Shelf
%     'wilkins' Wilkins Ice Shelf
%     'all' equivalent to bedmap2('shelves')
%
% outlineashelf(...,'LineProperty',LinePropertyValue) sets outline
% properties. 
%
% h = outlineashelf(...) returns the handle h of the shelf outline. 
%
% [h,lat,lon] = outlineashelf(...)
% 
%% Example 1
% Here we plot all ice shelfes and highlight Ross by making it red and plot
% FRIS in magenta. We'll also make Amery a green patch object. 
% 
% bedmap2('patchgl');
% outlineashelf('all'); % (equivalent to the command bedmap2('shelves'); )
% 
% [~,alat,alon] = outlineashelf('amery'); 
% patchm(alat,alon,'green') 
% scarlabel('Amery Ice Shelf');
% 
% outlineashelf('ross','color','red');
% scarlabel('Ross Ice Shelf','color','red'); 
% 
% outlineashelf('filchner','color','m','linewidth',1);
% scarlabel('Filchner Ice Shelf','color','m','fontweight','bold','rotation',30); 
% scarlabel('South Pole','horizontalalignment','left','marker','k*');
%
%% Example 2
% 
% figure
% bedmap2 'patchgl'
% mapzoom(-72,-85,'mapwidth',1500)
% 
% outlineashelf wilkins
% scarlabel('Wilkins Ice Shelf','color','blue')
% 
% outlineashelf('Abbot','color','g','linewidth',3)
% scarlabel('Abbot Ice Shelf','color',[.3 .5 .2],...
%     'fontweight','bold')
% 
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
% Created by Chad A. Greene
% The University of Texas at Austin
% August 2013.  Slightly updated October 2016. 
% 
% See also scarloc, scarlabel, and bedmap2. 

assert(exist('bedmap2shelf.mat','file')==2,'Cannot find the correct Bedmap2 data.') 

load bedmap2shelf.mat
plotall = false; % by default, don't plot all ice shelves. 
h = 1; 

switch lower(shelf)
    case {'filchner-ronne','fris','filchner','filchner ice shelf','ronne',...
            'ronne ice shelf'}; n = [1;16;21;35;46]; 
        disp('Plotting Filchner-Ronne Ice Shelf together'); 
    case {'ross','ross ice shelf'}; n = [2;23]; 
    case {'larsen c','larsen c ice shelf','larsen d','larsen d ice shelf'}; n = [3;48]; 
        disp('Plotting Larsen C and Larsen D together'); 
    case {'baudouin','baudouin ice shelf','borchgrevink','borchgrevink ice shelf',...
            'lazarev','lazarev ice shelf'}; n = 4; 
        disp('Plotting Baudouin, Borchgrevink, and Lazarev Ice Shelves together.'); 
    case {'fimbul','fimbul ice shelf','fimbulisen','jelbart','jelbart ice shelf','atka',...
            'atka ice shelf','ekstrom','ekstrom ice shelf'}; n = [5;36]; 
        disp('Plotting Fimbul, Jelbart, Atka, and Ekstrom ice shelves together'); 
    case {'riiser-larsen','riiser-larsen ice shelf','riiser','stancomb ice shelf',...
            'stancomb','stancomb brunt','stancomb-brunt','stancomb brunt ice shelf',...
            'stancomb-brunt ice shelf','brunt','brunt ice shelf'}; n = 6; 
        disp('Plotting Stancomb Brunt and Riiser-Larsen ice shelves together'); 
    case {'amery','amery ice shelf'}; n = 7;
    case {'getz','getz ice shelf'}; n = 8; 
    case {'abbot','abbot ice shelf'}; n = 9; 
    case {'george vi','george vi ice shelf'}; n = 10; 
    case {'sulzberger','sulzberger ice shelf'}; n = 11;
    case {'shackleton ice shelf','shackleton'}; n = 12; 
    case {'tracy','tremenchus','tracy ice shelf','tremenchus ice shelf'}; n = 12;
        disp('Plotting Tracy/Tremenchus with Shackleton'); 
    case {'west ice shelf','west'}; n = 13;
    case {'totten','totten glacier','totten glacier tongue'}; n = 14;
    case {'wilkins','wilkins ice shelf'}; n = 15; 
    case {'nickerson','nickerson ice shelf'}; n = 17;
        disp('Also plotting the adjacent Land Ice Shelf.')
    case {'land','land ice shelf'}; n = 17;
        disp('Also plotting the adjacent Nickerson Ice Shelf.')
    case {'stange','stange ice shelf'}; n = 18; 
    case {'dotson','dotson ice shelf'}; n = 19
        disp('Also plotting the adjacent Crosson Ice Shelf.')
    case {'crosson','crosson ice shelf'}; n = 19;        
        disp('Also plotting the adjacent Dotson Ice Shelf.')
    case {'prince harald','prince harald ice shelf'}; n = 20; 
    case {'pine island','pig','pine island glacier','pine island ice shelf',...
            'pine island bay','pine island bay ice shelf'}; n = 22;
    case {'bach','bach ice shelf','ill be bach'}; n = 24; 
    case {'mertz','mertz glacier tongue'}; n = 25;
    case {'moscow','muis','moscow university','moscow university ice shelf'}; n = 26;
    case {'drygalski','drygalksi ice tongue','drygalksi ice shelf'}; n = 27;
    case {'nivl','nivl ice shelf','nivlisen'}; n = 28; 
    case {'thwaites','thwaites glacier tongue'}; n = 29;
    case {'mariner','mariner ice shelf','mariner glacier tongue'}; n = 30;     
    case {'rennick','rennick ice shelf','rennick bay ice shelf'}; n = 31;
    case {'ninnis','ninnis glacier tongue'}; n = 32;
    case {'cook','cook ice shelf'}; n = 33;
    case {'nansen','nansen ice shelf'}; n = 34; 
    case {'cosgvove','cosgrove ice shelf'}; n = 37; 
    case {'venable','venable ice shelf'}; n = 38; 
    case {'larsen b','larsen b ice shelf'}; n = 39;
    case {'publications','publications ice shelf'}; n = 40; 
    case {'vigrid','vigridisen','vigrid ice shelf'}; n = 41; 
    case {'holmes','holmes glacier'}; n = 42;
    case {'larsen e','larsen e ice shelf'}; n = 43; 
    case {'swinburne','swinburne ice shelf'}; n = 44; 
    case {'conger','conger glacier','conger ice shelf'}; n = 45;
    case {'vincennes','vincennes bay'}; n = 47;
    case {'quar','quarisen','quar ice shelf'}; n = 49; 
    case {'rayner','rayner ice shelf','rayner glacier','thyer','thyer glacier',...
            'thyer ice shelf','rayner-thyer','rayner thyer'}; n = 50; 
    case {'dibble','dibble glacier','dibble glacier tongue'}; n = 56;
    case {'all'}; plotall=true; 

end

if ~strcmpi(shelf,'all'); 
lat = shlflat{min(n)};
lon = shlflon{min(n)};
end
          
          
if license('test','map_toolbox')
   if plotall
       try
       h = NaN(604,1); 
       for n = 1:604
           h(n) = plotm(shlflat{n},shlflon{n});
       end
       end
   else 
       try
       h = NaN(size(n));
       for nn = 1:length(n)
           h(nn) = plotm(shlflat{n(nn)},shlflon{n(nn)});
       end
       end
       lat = shlflat{n(min(nn))};
       lon = shlflon{n(min(nn))};
   end

   % Set varargin inputs: 
   for k = 1:2:length(varargin)
       set(h,varargin{k},varargin{k+1}); 
   end

end

% delete the outputs if they are not requested: 
if nargout==0
    clear h lat lon; 
end

end

