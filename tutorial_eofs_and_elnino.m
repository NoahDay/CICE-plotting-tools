close all
load pacific_sst.mat

%whos

datestr(t([1,end]));

mean(diff(t));


[Lon,Lat] = meshgrid(lon,lat);

% El nino southern oscilation index
idx = enso(sst,t,Lat,Lon);

anomaly(t,idx)



%% Plot power density spectrum
load train
t = (0:length(y)-1)/Fs;

plot(t,y)
box off
xlabel 'time (s)'
plotpsd(y,Fs)
xlabel 'frequency (Hz)'

plotpsd(idx,12) 
plotpsd(single(idx),12,'logx','lambda')
xlabel('Periodicity (years)')

%%
% Calculate the mean SST
sst_mean = mean(sst,3);
imagescn(lon,lat,sst_mean)
cb = colorbar;

% Calculate trend (deg/yr)
sst_trend = 365.25*trend(sst,t,3);
imagescn(lon,lat,10*sst_trend)
cb = colorbar;
ylabel(cb,'temperature trend {\circ}C per decade')
cmocean('balance','pivot')

%% Eofs identify the modes of variability of a system

imagescn(lon,lat,eof(sst,1))
colorbar;
cmocean('balance','pivot')
title 'eof first mode'

%% Mark a location of interest
hold on
plot(lon(12),lat(10),'ks')
hold off
%%
% Get time series at that location of interest
sst1 = squeeze(sst(10,12,:));

plot(t,sst1)
datetick
% Deaseason the data
sst1_ds = deseason(sst1,t);

hold on
plot(t,sst1_ds)
hold off
%%
sst_ds = deseason(sst,t);

imagescn(lon,lat,eof(sst_ds,1))
colorbar
cmocean('balance','pivot')

%%
sst_ds_dt = detrend3(sst_ds);

sst_anom_var = var(sst_ds_dt,[],3); % alont the third dimension

imagescn(lon,lat,sst_anom_var)
caxis([0,1])
%% What are the modes of variability in this deseasoned detrended data

% Calculate eofs
[eof_maps,pc,expv] = eof(sst_ds_dt,6);
% PC is the principal component time series
clf
subplot(3,2,1)
imagescn(lon,lat,eof_maps(:,:,1))
axis off
cmocean('balance','pivot')
axis image

subplot(3,2,2)
plot(t,pc(1,:))
axis tight
box off
datetick

subplot(3,2,3)
imagescn(lon,lat,eof_maps(:,:,2))
axis off
cmocean('balance','pivot')
axis image

subplot(3,2,4)
plot(t,pc(2,:))
axis tight
box off
datetick

subplot(3,2,5)
imagescn(lon,lat,eof_maps(:,:,3))
axis off
cmocean('balance','pivot')
axis image

subplot(3,2,6)
plot(t,pc(3,:))
axis tight
box off
datetick

sgtitle 'The first three principal components'
%%
expv
%%
clf
subplot(1,2,1)
h1 = imagescn(lon,lat,sst_ds_dt(:,:,1));
title 'observed sst anomaly'
cmocean bal
caxis([-1,1]*2.5)

subplot(1,2,2)
h1 = imagescn(lon,lat,eof_maps(:,:,1).*pc(1,1));
title 'reconstructed sst anomaly'
cmocean bal
caxis([-1,1]*2.5)

sgtitle(datestr(t(1)))
%%

h2.CData = eof_maps(:,:,1)*pc(1,1) + ...
           eof_maps(:,:,2)*pc(2,1) + ...
           eof_maps(:,:,3)*pc(3,1) + ...
           eof_maps(:,:,4)*pc(4,1) + ...
           eof_maps(:,:,5)*pc(5,1) + ...
           eof_maps(:,:,6)*pc(6,1);
% Reconstruct sst anomalies from first 5 modes
sst_ds_dt_r = reof(eof_maps,pc,1:5);
for k = 1:120
    h1.CData = sst_ds_dt(:,:,k);
    h2.CData = sst_ds_dt_r(:,:,k);

    subplot(1,2,1)
    h1 = imagescn(lon,lat,h1.CData );
    title 'observed sst anomaly'
    cmocean bal
    caxis([-1,1]*2.5)
    
    subplot(1,2,2)
    h1 = imagescn(lon,lat,h2.CData);
    title 'reconstructed sst anomaly'
    cmocean bal
    caxis([-1,1]*2.5)
    pause(0.1)
    sgtitle(datestr(t(k)))
end

