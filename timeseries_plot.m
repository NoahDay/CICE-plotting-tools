%% Create time series from run files
clear
close
% fid = fopen("cases/3month/logs/cice.runlog.220127-122445")
%  Data = fread(fid);
%  CharData = char(Data);
%  while 1
%             tline = fgetl(fid);
%             if ~ischar(tline), break, end
%             disp(tline)
% end
%  
%  fclose(fid);
%  %disp(Data);

ice_extent_3month = [9352691.727138786, 10417525.631410243, 11552288.354975263, 12438142.707410125, 13146347.793605084, 13591250.059821103, 13955683.031831602, 14303909.124725379, 14577858.748199204, 14959615.946081417, 15306485.061280703, 15695954.35552291, 16008562.091676988, 16215486.025221244, 16491932.056667447, 16791689.373903196, 17034680.54856902, 17283164.91086511, 17484331.25655016, 17817019.2700459, 18165826.846462812, 18510506.15149546, 18792757.368977223, 18974432.439918935, 19172523.08904214, 19391903.98906456, 19531649.302846745, 19675797.312850136, 19783329.93242231, 19937675.293143447, 20216670.69061112, 20467200.18920361, 20766509.026904676, 21036983.713430796, 21282016.918412652, 21582023.25523878, 21743316.499327645, 21804835.990868747, 22028363.293441903, 22122272.862666897, 22051063.783808723, 22092318.68040637, 22207603.482449796, 22481691.500704408, 22724767.60705103, 22866975.623221364, 22967468.07140107, 23010790.512440026, 23067534.292171687, 23192565.19201289, 23276172.03893974, 23238317.559958726, 23208257.952295974, 23234979.823786404, 23303239.66052915, 23266649.506475512, 23321347.03114365, 23317047.06435271, 23340605.065086797, 23280631.734995075, 23281984.542034913, 23240908.76619559, 23284319.222907756, 23226772.779761985, 23265998.391622208, 23270308.985865988, 23242162.60150222, 23320438.09672123, 23395318.990520652, 23403083.833157536, 23436184.650607355, 23378643.43887168, 23428103.375837203, 23435837.184564404, 23537789.21218092, 23556257.345424537, 23415123.439880617, 23416581.023139156, 23436377.819403987, 23408491.907128245, 23351418.616395842, 23262022.77773243, 23176228.913401015, 23110362.815987825, 23097712.84900041, 23143921.624585334, 23095399.28154743, 22964047.05800472, 22935125.520272758, 23050948.64665129, 23055118.002206624, 23030868.40298112];
ice_extent_1year = [9117580.534278475, 8168877.212313742, 7360205.811863311, 6757331.843435727, 6283689.745810203, 5862781.760577699, 5474957.532287999, 5082160.915263696, 4775879.630731631, 4520538.11310633, 4279995.869313685, 4075412.9033696344, 3850297.653657133, 3681154.619478512, 3461135.032023103, 3307370.9167422825, 3135573.811110725, 2971648.36006258, 2813630.0045050457, 2682488.2866993, 2553929.8350452227, 2427520.9430597126, 2298519.26532907, 2152835.852986114, 2059421.5718364366, 1998669.4829396405, 1946743.1564285005, 1913943.9100719227, 1862752.5411649786, 1791251.1907794387, 1785332.7284562902, 1748501.3926486599, 1699977.0838345564, 1664317.5549368656, 1629561.9045807463, 1641308.223588646, 1624628.5258779323, 1623077.2223425547, 1610772.6901391428, 1625260.2992042664, 1610283.2057999515, 1617991.879827451, 1622576.9575092734, 1624735.6058663088, 1643638.0036419104, 1644502.4792805365, 1626848.401584588, 1660479.4749449277, 1675623.1287770367, 1706892.1889267531, 1738416.90835496, 1748077.6616184588, 1766982.936384959, 1759205.7968347566, 1771420.3819883543, 1803705.3440879444, 1832475.1374506573, 1877273.4605763084, 1921559.873701711, 2030320.0088560127, 2079735.1787323514, 2171196.4042702373, 2289830.750196537, 2345797.5457769, 2425844.225131105, 2523286.6936102547, 2575676.5164685166, 2647200.0915035866, 2745091.6551369857, 2831334.9961279677, 2851270.6856709863, 2925601.6683975565, 2959192.7183812773, 3024650.254407117, 3086747.40675553, 3187431.8162430944, 3303768.6636379063, 3415821.7662520073, 3516577.3967331788, 3609832.7304188213, 3708573.100082026, 3865452.4289497384, 3998315.36393668, 4170807.2627647333, 4332540.555328448, 4547410.709781474, 4709434.396935783, 4828070.543818521, 4940722.225691992, 5136807.049866636, 5419634.2348409165, 5614347.408761663, 5779752.040993309, 6048043.952973603, 6266030.450900743, 6389077.65206702, 6539110.684480219, 6677036.476066584, 6854316.233302162, 7066777.897248114, 7239571.964537569, 7353702.060212138, 7498976.898714985, 7645520.3488337295, 7807175.89485616, 7888030.540741754, 7950575.396093543, 8044878.632923863, 8231222.747256357, 8465999.591405218, 8582728.153624231, 8750322.637806488, 8889367.909313457, 8983700.89166821, 9141085.588315055, 9313317.753655884, 9590637.793934979, 9761588.799853746, 9937197.780031106, 10054441.671549298, 10167776.875448871, 10235265.955287017, 10303805.053772142, 10568536.03486194, 10805306.568160132, 11024734.509123605, 11262521.513951277, 11487195.551563803, 11772456.30996805, 11972092.726886012, 12174581.827440383, 12221155.26738832, 12419973.844417773, 12612754.256488316, 12742240.19316875, 12915286.279246476, 13051192.33845084, 13349539.871894008, 13593478.469158765, 13698649.315127784, 13846121.26323189, 14108554.26242919, 14456910.409018928, 14773255.653101597, 15153159.91036387, 15392512.806250697, 15571972.520875398, 15649200.915959746, 15698166.348184234, 15801254.27067525, 15903241.68564115, 16030866.621728983, 16153535.598257506, 16341193.642020708, 16449072.766360108, 16544210.411868257, 16742742.239186665, 16823843.003176764, 17015917.75981427, 17206713.228305444, 17439880.96192583, 17707572.87894086, 17945471.774613105, 18128169.59360975, 18361183.106503136, 18496687.826482587, 18697476.14703142, 18882408.107163686, 19258745.388909336, 19433156.32653561, 19540672.571435805, 19670888.61691527, 19793037.96931619, 20014392.553787522, 20199327.018928397, 20351081.207213156, 20502551.94029201, 20632242.767002076, 20813914.915351678, 20923758.04906957, 20963530.485764142, 21260407.635988805, 21487698.925089374, 21743233.688824333, 21765693.65436347, 21923740.644057646, 22062357.827975363, 22172775.324893955, 22339731.940096088, 22492616.72951843, 22669340.50782111, 22836188.46398322, 23049699.030415546, 23175067.280440662, 23271695.11240553, 23275404.9310422, 23436178.316828545, 23506600.095064607, 23528981.42168121, 23635297.51110299, 23792085.548064098, 23916042.074410688, 23950181.97220959, 24057340.77866728, 24078565.374860987, 24225495.01153596, 24374213.149928153, 24440584.43920986, 24488210.59573439, 24584502.842427462, 24580128.02448515, 24705759.83686371, 24906259.61124762, 25046994.780327905, 25267100.80466559, 25415820.886427283, 25479313.0740237, 25515665.747519303, 25505372.441396043, 25603191.424969234, 25559483.418592416, 25422657.7051659, 25303956.700919773, 25294472.827103894, 25490679.118370023, 25635014.050358374, 25686671.020465896, 25743448.79367877, 25707801.983533476, 25712137.874409612, 25781317.913290188, 25685183.57805953, 25594053.817731146, 25525311.641522516, 25405869.27405336, 25400188.89246085, 25401053.057659462, 25366996.67364332, 25313020.547098614, 25234923.220764406, 25221510.910063457, 25151445.50546411, 25019354.018299807, 24938585.690650046, 24882615.42071793, 24878732.36515478, 24836889.7893421, 24815259.431372996, 24872883.98767664, 24935752.511847787, 24873842.523713477, 24901932.49069651, 24850957.262335185, 24846984.899863508, 24847735.11565934, 24901136.851003256, 24898038.834694553, 24724025.557993647, 24705093.55068495, 24648477.356786966, 24593591.295971833, 24562186.242144227, 24399532.460118927, 24294160.399402, 24252078.888654616, 24172443.79998293, 24139330.541231137, 24052761.54439807, 23918426.833527286, 23895537.371639665, 23959675.343262266, 23958054.019735746, 23998185.977822825, 24032958.67665739, 23963657.95346431, 23847572.08499613, 23635246.94428325, 23541586.454655558, 23451530.161745347, 23358113.14971795, 23185928.65007764, 23010787.080760077, 22861839.512695797, 22796035.55255371, 22676639.19394987, 22625281.48970038, 22538814.42193084, 22415209.782925807, 22257451.848334387, 22004715.60879959, 21886612.91275928, 21763469.774490874, 21702796.139569722, 21572330.90293304, 21506493.072596867, 21366868.776025347, 21296402.33154557, 21269312.421571027, 21099761.825065188, 20999981.90870933, 20924976.980821695, 20825649.228824243, 20659013.97625745, 20517701.221166465, 20382223.00501972, 20278654.0829101, 20168527.104154162, 20062604.607149895, 19820724.509366084, 19684504.182176545, 19515379.23997857, 19333725.850399576, 19154114.007701255, 19040172.514539976, 18772904.715563882, 18685150.72522471, 18565676.174327414, 18336332.84847453, 18174810.198329534, 18058079.686783604, 17940850.988295328, 17832929.086838186, 17718290.052026078, 17613681.54480148, 17500981.75769098, 17405566.606891423, 17325819.53168255, 17187510.6186702, 17068100.949448656, 16933805.551647674, 16823992.839604124, 16712048.582391951, 16593951.52354773, 16502572.472548628, 16435337.643845726, 16312978.61358774, 16222024.519079847, 16130073.297620656, 15932443.696707387, 15802486.223549927, 15698321.189254412, 15544203.0840512, 15472861.688928375, 15358608.209294464, 15252672.038146254, 15137062.444274167, 15039142.509585941, 14895423.760290155, 14779053.849100715, 14674427.52758039, 14571451.216188012, 14469791.470212378, 14371423.900665147, 14231651.266580043, 14091727.174621956, 13963813.503104242, 13762120.087932602, 13641807.402106369, 13499991.472435908, 13330608.02289916, 13162789.021724738, 12966485.8164872, 12822090.390381722, 12628718.915873105, 12339885.666513892];
ice_extent_3month_july  = [14815667.0405731, 15051859.309609888, 15296744.344786782, 15562524.244159242, 15879316.548869094, 16235754.823892826, 16587810.27182586, 16889681.45169224, 17064205.064946566, 17381213.46902154, 17733302.07932438, 18119284.414009493, 18441172.864103742, 18623936.480156414, 18992962.538384285, 19352077.217574995, 19524729.79413071, 19683068.50244795, 19900532.519975003, 20151574.803932503, 20391437.171559315, 20680016.805355188, 20897881.55341655, 21038227.480033908, 21198806.03738401, 21406265.31014319, 21549715.689209577, 21697484.07743515, 21794134.220510405, 21899047.700126275, 22176158.223762173, 22378714.997968458, 22689096.817698848, 22896798.212792236, 23094959.489657205, 23370730.042232793, 23475985.834896058, 23500219.21324057, 23690496.203192536, 23690962.951884285, 23675210.000377323, 23570657.340542316, 23652337.41975308, 23821434.40174814, 24024433.338724274, 24081967.065249737, 24155357.75513813, 24237969.888017986, 24191416.71374074, 24317586.024576847, 24344919.676202543, 24334709.1592899, 24209444.84035634, 24142618.665842045, 24185428.05485896, 24214337.93112205, 24217817.797794923, 24177090.069624763, 24124803.920616273, 24154852.98146663, 24090146.365663394, 24012650.33747853, 23987388.535630554, 23981285.17086173, 23922974.55293774, 23938166.577633817, 23854337.386226807, 23963872.913599823, 23990462.65797366, 23978286.332591757, 23971506.136661205, 23981234.065616705, 23988693.965513337, 24034073.5230438, 24112107.09110873, 24078924.180460136, 23972162.876964234, 23964198.827436954, 23971282.365792237, 23918156.392722666, 23827235.782582577, 23722343.680749208, 23603058.16316951, 23560094.59743566, 23529514.822263964, 23521296.985201914, 23441811.4992141, 23293641.448797207, 23294840.357052192, 23421875.706258222, 23452102.616009865, 23431318.119225964];
ice_extent_6month = [4240193.956579872, 4578493.649003474, 4891051.431613879, 5300671.914499691, 5575393.246823255, 5782074.402084323, 5951339.052709573, 6108656.07223593, 6334879.501265983, 6578121.46855489, 6812287.761379639, 6975193.101283698, 7200576.05652514, 7397680.683564755, 7541317.687652813, 7682714.929279335, 7769268.203017274, 7838322.657678657, 8068142.028874659, 8273660.624742303, 8439249.73446856, 8612353.943967171, 8703745.775317458, 8786286.046511993, 8998499.738147493, 9187133.997548377, 9475087.410195814, 9638442.126869595, 9839122.045432974, 10003436.616010485, 10142180.90009644, 10236589.534816459, 10302538.5297649, 10471619.933776263, 10714911.680389335, 10950313.7903804, 11243332.961990433, 11519483.348821022, 11732845.570050571, 11931956.302420748, 12095942.856788816, 12215675.520023577, 12418943.400933126, 12639438.850832758, 12760226.10006395, 12925529.081240593, 13155313.05386244, 13416948.171004264, 13692107.82292878, 13832005.328336539, 13954721.107664686, 14151317.20962447, 14509447.031262303, 14892505.795076903, 15202233.92249457, 15416105.864967965, 15575547.91309188, 15701951.173271023, 15709004.12711046, 15798255.18073421];
%% 3 month jan initial conditions
t1 = datetime(2005,7,1);
t2 = datetime(2005,9,30);
dates_3m = datevec(t1:t2);

ts1 = timeseries(ice_extent_3month,1:length(ice_extent_3month));
ts1.Name = 'Antarctic ice extent (km^2)';
ts1.TimeInfo.Units = 'days';
ts1.TimeInfo.StartDate = '01-Jul-2005';     % Set start date.
ts1.TimeInfo.Format = 'mmm dd, yy';       % Set format for display on x-axis.


%% 1 year
t1 = datetime(2005,1,1);
t2 = datetime(2005,12,31);
dates_1y = datevec(t1:t2);

ts2 = timeseries(ice_extent_1year,1:length(ice_extent_1year));
ts2.Name = 'Antarctic ice extent (km^2)';
ts2.TimeInfo.Units = 'days';
ts2.TimeInfo.StartDate = '01-Jan-2005';     % Set start date.
ts2.TimeInfo.Format = 'mmm dd, yy';       % Set format for display on x-axis.

ts2.Time = ts2.Time - ts2.Time(1);        % Express time relative to the start date.

%% 3 month july init
t1 = datetime(2005,7,1);
t2 = datetime(2005,9,30);
dates_1y = datevec(t1:t2);

ts3 = timeseries(ice_extent_3month_july,1:length(ice_extent_3month_july));
ts3.Name = 'Antarctic ice extent (km^2)';
ts3.TimeInfo.Units = 'days';
ts3.TimeInfo.StartDate = '01-Jul-2005';     % Set start date.
ts3.TimeInfo.Format = 'mmm dd, yy';       % Set format for display on x-axis.

ts3.Time = ts3.Time - ts3.Time(1);        % Express time relative to the start date.

%% 6 month
t1 = datetime(2005,4,1);
t2 = datetime(2005,5,31);
dates_1y = datevec(t1:t2);

ts4 = timeseries(ice_extent_6month,1:length(ice_extent_6month));
ts4.Name = 'Antarctic ice extent';
ts4.TimeInfo.Units = 'days';
ts4.TimeInfo.StartDate = '01-Apr-2005';     % Set start date.
ts4.TimeInfo.Format = 'mmm dd, yy';       % Set format for display on x-axis.

ts4.Time = ts4.Time - ts4.Time(1);        % Express time relative to the start date.

%% Plotting
line_width = 3;
plot(ts1,'LineWidth',line_width)
hold on
plot(ts2,'LineWidth',line_width)
plot(ts3,'LineWidth',line_width)
plot(ts4,'LineWidth',line_width)
    lgd=legend('Winter 1','1 year','Winter 2', 'AprilMay','Position',[0.21 0.7 0.1 0.2]);
    lgd.FontSize = 14;
hold off