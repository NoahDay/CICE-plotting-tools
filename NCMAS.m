%runtime_per_step = [0.42, 0.19, 0.11, 0.067];
data = [98.95131202855121, 0.42569481518729646;
        181.45392605755094, 0.22835291367858995;
        346.0199742669309, 0.11740109055380506;
        647.0546695999008, 0.06297676208532169;
        1097.2505626955144, 0.042973103690445213;
        1934.9094096180622, 0.031031536417971687];

% 96 (2), 192 (4), 384 (8), 768 (16), 1152 (24), 2304 (48)
%runtime_per_step = [runtime_per_step.*1.19*1.026];
num_cpu = [100, 210, 400, 800];
Fit = polyfit(log10(data(1:4,1)),log10(data(1:4,2)),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 

close all
addpath functions
conFigure(11)
f = figure;
scatter(data(:,1),data(:,2),'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
hold on
ord = polyval(Fit,log10(data(:,1)));
fit = 10.^ord;
plot(data(:,1),fit,'--')

set(gca,'xscale','log')
set(gca,'yscale','log')

xlim([50,0.2*10^4])
ylim([10^(-2),1])
xlabel('Number of CPUs')
ylabel('Runtime/step [s]')
xticks([10,100,200,400,1000,2000]);
xticklabels({'10','100','200','400','1000','2000'})




exportgraphics(f,'runtime.pdf','ContentType','vector')

%%

data = [2, 0.42569481518729646;
        4, 0.22835291367858995;
        8, 0.11740109055380506;
        16, 0.06297676208532169;
        24, 0.048973103690445213;
        48, 0.031031536417971687];


num_cpu = [100, 210, 400, 800];
Fit = polyfit(log10(data(:,1)),log10(data(:,2)),1); % x = x data, y = y data, 1 = order of the polynomial i.e a straight line 
Fit(1) = -1;

Fit(2) = -0.046;
close all
addpath functions
conFigure(11)
f = figure;
scatter(data(:,1),data(:,2),80,'MarkerEdgeColor',[0 .5 .5],...
              'MarkerFaceColor',[0 .7 .7],...
              'LineWidth',1.5)
hold on
ord = polyval(Fit,log10(data(:,1)));
fit = 10.^ord;
plot(data(:,1),fit,'--')

set(gca,'xscale','log')
set(gca,'yscale','log')

%xlim([50,0.2*10^4])
ylim([10^(-2),1])
xlabel('Number of nodes (48 cores each)')
ylabel('Runtime/step [s]')
xticks([2,4,8,16,24,48]);
xticklabels({'2','4','8','16','24','48'})




exportgraphics(f,'runtime_node.pdf','ContentType','vector')