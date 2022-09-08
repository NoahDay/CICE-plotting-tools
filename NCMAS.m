runtime_per_step = [0.42, 0.19, 0.11, 0.067];
runtime_per_step = runtime_per_step.*1.19*1.026;
num_cpu = [100, 210, 400, 800];
Fit = polyfit(log10(num_cpu),log10(runtime_per_step),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 

close all
addpath functions
conFigure(11)
f = figure;
scatter(num_cpu,runtime_per_step,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
hold on
ord = polyval(Fit,log10(num_cpu));
fit = 10.^ord;
plot(num_cpu,fit,'--')

set(gca,'xscale','log')
set(gca,'yscale','log')

xlim([10,0.2*10^4])
ylim([10^(-2),1])
xlabel('Number of CPUs')
ylabel('Runtime/step [s]')
xticks([10,100,200,400,1000,2000]);
xticklabels({'10','100','200','400','1000','2000'})




exportgraphics(f,'runtime.pdf','ContentType','vector')