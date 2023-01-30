clear all
close all
clc
addpath functions

SIC = 0.15; 
filename = "/Users/a1724548/github/cice-dirs/runs/initial/history/iceh_inst.2005-01-01-03600.nc";
filename = "/Volumes/NoahDay5TB/WIMonAlessandroRun/history/iceh.2017-09-01.nc";
% Read the header
ncdisp(filename)
% 
grid = 'om2';

nw = 16;
nfreq = nw;
fmin = 1/1000;%1/50;%1/16; % freq min
fmax = 1;%1/2;%1/6; % freq max
om1=2*pi*fmin; % ang freqs, rad/s
om2=2*pi*fmax;
om_0 = (om2 - om1)/(nw-1); % steps
gravit = 9.81;
% Calculate wave numbers and wavelengths
for lp_i=1:nw
  omega(lp_i)        = om1 + (lp_i-1)*om_0; % Frequency, rad/s
  T(lp_i)            = 2*pi/omega(lp_i); % Period, s
  lam_wtr_in(lp_i)   = gravit*(T(lp_i)^2)/2/pi; % Wave number
  k_wtr_in(lp_i)     = 2*pi/lam_wtr_in(lp_i); % Wave length
end


% Importing om values
for lp_i=1:nfreq
    om(lp_i)        = om1 + (lp_i-1)*om_0;
end

dwavefreq(1) = om(1);
for lp_i=2:nfreq
  dwavefreq(lp_i) = om(lp_i) - om(lp_i-1);
end 


%cd /Users/a1724548/Github/CICE-plotting-tools
%wave_spectrum = data_format(filename,"wave_spectrum");
%for i = 1:321
%    for j = 1:384
%        twodspec(i,j) = max(wave_spectrum(i,j,:))>eps;
%    end
%end

dfreq = dwavefreq/(2*pi); % rad/s to 1/s
for i = 1:321
    for j = 1:384
        work(:) = wave_spectrum(i,j,:);%/(2*pi); % m^2/s
        hs(i,j) = sqrt(work*dwavefreq'); % sqrt(m^2/s*(1/s)^2)
        smom0 = simpson_integrate(wave_spectrum(i,j,:)/(2*pi), 0, om, nw);
        hm0(i,j)= 4*sum(sqrt(smom0));
    end
end
f = figure;
[p,a] = map_plot(hs,"wave_sig_ht","SH","gx1",[0,max(max(hs))]);
a.Label.String = "$H_s$ (m)";
a.Label.Interpreter = "latex";
title('Reimann integral')
 %exportgraphics(f3,'fsdrad2.pdf','ContentType','vector')
colormap parula

wave_sig_ht = data_format(filename,"wave_sig_ht");
f = figure;
[p,a] = map_plot(wave_sig_ht,"wave_sig_ht","SH","gx1",[0,max(max(wave_sig_ht))]);
a.Label.String = "$H_s$ (m)";
a.Label.Interpreter = "latex";
title('CICE output')
 %exportgraphics(f3,'fsdrad2.pdf','ContentType','vector')
colormap parula

f = figure;
hm0 = hm0/(2*pi);
[p,a] = map_plot(hm0,"wave_sig_ht","SH","gx1",[0,max(max(hm0))]);
a.Label.String = "$H_s$ (m)";
a.Label.Interpreter = "latex";
title("Simpson's integral")
 %exportgraphics(f3,'fsdrad2.pdf','ContentType','vector')
colormap parula

ran = 2;
f = figure;
[p,a] = map_plot(hs-wave_sig_ht,"wave_sig_ht","SH","gx1",[-ran,ran]);
a.Label.String = "$H_s$ difference (m)";
a.Label.Interpreter = "latex";
title("Difference ($\sum S(f) \Delta S(f)$ - CICE output)",'Interpreter','latex')
 %exportgraphics(f3,'fsdrad2.pdf','ContentType','vector')
cmocean('balance',11)

f = figure;
[p,a] = map_plot(hs-hm0,"wave_sig_ht","SH","gx1",[-ran,ran]);
a.Label.String = "$H_s$ difference (m)";
a.Label.Interpreter = "latex";
title("Difference ($\sum S(f) \Delta S(f) - \texttt{simp}$)",'Interpreter','latex')
 %exportgraphics(f3,'fsdrad2.pdf','ContentType','vector')
cmocean('balance',11)

f = figure;
[p,a] = map_plot(hm0-wave_sig_ht,"wave_sig_ht","SH","gx1",[-ran,ran]);
a.Label.String = "$H_s$ difference (m)";
a.Label.Interpreter = "latex";
title("Difference ($\texttt{simp}$ - CICE output)",'Interpreter','latex')
 %exportgraphics(f3,'fsdrad2.pdf','ContentType','vector')
cmocean('balance',11)

%% Spectral moment
clc
S(:) = wave_spectrum(1,30,:);
wave_sig_ht(1,30)
dum_S = S;
mom = 0;
om_in = om;
nw_in = nw;


domega = omega - [0,omega(1:15)];
(dum_S*2*pi)*(domega*2*pi)'

smom0 = simpson_integrate(dum_S, mom, om_in, nw_in);

4*sum(sqrt(smom0))

function dum_sm0 = simpson_integrate(dum_S, mom, om_in, nw_in)
    dum_simp(1) = 1;
    dum_simp(nw_in) = 1;
    nth_in = 1;
    dom_local = om_in(2)-om_in(1);
     if nth_in == 1
      dth_local = 3;         % for Simpson's rule
     else
      dth_local = th_in(2)-th_in(1);
     end

    for loop_w=2:2:nw_in-1
      dum_simp(loop_w) = 4;
     end 

    for loop_w=3:2:nw_in-1
      dum_simp(loop_w) = 2;
     end 

     dum_simp_th(1) = 1;
     dum_simp_th(nth_in) = 1;

     for loop_th=2:2:nth_in-1
      dum_simp_th(loop_th) = 4;
     end

    for loop_th=3:2:nth_in-1
      dum_simp_th(loop_th) = 2;
    end

    for loop_w=1:nw_in
      for loop_th=1:nth_in
       wt_simp(loop_w+nw_in*(loop_th-1)) = dum_simp(loop_w)*dum_simp_th(loop_th);
      end 
     end 

     for loop_w=1:nw_in
      for loop_th=1:nth_in
       wt_int(loop_w+nw_in*(loop_th-1)) = (dom_local/3d0)*(dth_local/3d0)*wt_simp(loop_w+nw_in*(loop_th-1));
       dum_v(loop_w+nw_in*(loop_th-1))  = dum_S(loop_w+nw_in*(loop_th-1))*(om_in(loop_w)^mom);
      end 
     end 

    dum_sm0   = wt_int.*dum_v;


end
