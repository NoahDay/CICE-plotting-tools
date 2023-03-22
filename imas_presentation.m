addpath functions
%clear all
% Truncated power law
% N(x|d_{min} <= x <= l_{max}) = Cx^{-alpha}
alpha = 2

%d_min = [1,100,500];
%l_max = [100,500,1000];
alpha = [1.5, 2, 3];
d_min = 2.5;
l_max = 10000;
n = 1000;
p_fsd = zeros(3,n);
for i = 1:3
    x = linspace(d_min,l_max,n);
    f =  @(x) x.^(-alpha(i));
    y = f(x);
    Int = cumtrapz(x,y);
    Intv = @(a,b) max(Int(x<=b)) - min(Int(x>=a));
    A = Intv(d_min, l_max);
    C = 1/A;
    p_fsd(i,:) = y.*C;
end

for i = 1:length(alpha)
   label_vec{i} = strcat('$\alpha$',sprintf(' = %g',alpha(i)));
end
label_vec = ["Consolidated","Unconsolidated","Ice edge"];
close all  
conFigure(30)
f = figure;
loglog(linspace(d_min,l_max,n),(p_fsd),'Linewidth',5)
yticks('auto')
ylabel('Probability')
xlabel('Floe radius [m]')
legend(label_vec,'Interpreter','latex')
title('Power law')
exportgraphics(f,'powerlaw.pdf','ContentType','vector')

%%
label_vec = ["$r_{\max} \to \infty$","$r_{\max} = 850$"];
[minValue,closestIndex] = min(abs(850-x));
i = 1;
int_extra = sum(p_fsd(i,closestIndex+1:end));
p_fsd_extra = p_fsd(1,1:closestIndex);
p_fsd_extra(end) = p_fsd_extra(end)+int_extra;


conFigure(30)
f = figure;
loglog(x,p_fsd(1,:),'Linewidth',5)
hold on
loglog(x(1:closestIndex),p_fsd_extra(1,:),'Linewidth',5)
hold off
yticks('auto')
ylabel('Probability')
xlabel('Floe radius [m]')
legend(label_vec,'Interpreter','latex')
exportgraphics(f,'fsd_extra.pdf','ContentType','vector')

%%
clf
alpha = 2;
d_min = 2.5;
l_max = 10;
n = 12;
p_fsd = zeros(1,n);
x = linspace(d_min,l_max,n);
f =  @(x) x.^(-alpha);
y = f(x);
Int = cumtrapz(x,y);
Intv = @(a,b) max(Int(x<=b)) - min(Int(x>=a));
A = Intv(d_min, l_max);
C = 1/A;
p_fsd(1,:) = y.*C;

%for i = 1:length(alpha)
%   label_vec{i} = strcat('$\alpha$',sprintf(' = %g',alpha(i)));
%end
label_vec = ["Consolidated","Unconsolidated","Ice edge"];
close all
conFigure(30)
f = figure;
%plot(linspace(d_min,l_max,n),(p_fsd),'Linewidth',5)
bar(1:n,(p_fsd))
yticks('manual')
xticks('manual')
ylabel('Probability')
xlabel('Floe category')
font_size = 25;
% One arrow from left to right with text on left side
x = [0.74 0.79];    % adjust length and location of arrow 
y = [0.8 0.8];      % adjust hieght and width of arrow
annotation('textarrow',x,y,'String',' Lateral growth ','FontSize',font_size,'Linewidth',2,'Interpreter','latex')
Xadj = 1.34;
x = [0.74 0.79];
y = [0.7 0.7];
annotation('textarrow',-x+Xadj,y,'String',' Lateral melt ','FontSize',font_size,'Linewidth',2,'Interpreter','latex')
% Arrow with two head at both end and text between
y = [0.6 0.6];    
Xadj = 1.34;      % adjust location of left arrow starting point (the sum of this with 'x' should not be negative)
annotation('textarrow',x,y,'String',' New floes  ','FontSize',font_size,'Linewidth',2,'Interpreter','latex')
annotation('textarrow',-x+Xadj,y,'String','','FontSize',14,'Linewidth',2)
% Welding
x = [0.74 0.79];
y = [0.5 0.5];
annotation('textarrow',x,y,'String',' Welding ','FontSize',font_size,'Linewidth',2,'Interpreter','latex')
% Wave breakup
Xadj = 1.34;
x = [0.74 0.79];
y = [0.4 0.4];
annotation('textarrow',-x+Xadj,y,'String',' Wave breakup ','FontSize',font_size,'Linewidth',2,'Interpreter','latex')


exportgraphics(f,'change_fsd.pdf','ContentType','vector')
%% joint FSTD

% Truncated power law
% N(x|d_{min} <= x <= l_{max}) = Cx^{-alpha}
alpha = 2

%d_min = [1,100,500];
%l_max = [100,500,1000];
alpha = linspace(3,-1,5);
alpha = [3, 2.5, 2, 0.5, -2];
d_min = 2.5;
l_max = 10000;
ncat = 5;
nfsd = 12;
p_fsd = zeros(ncat,n);
weights = [0.11, 0.14, 0.18, 0.02, 0.0001];
for i = 1:ncat
    %x = linspace(d_min,l_max,12);
    x = 1:12;
    f =  @(x) x.^(-alpha(i));
    y = f(x);
    p_fsd(i,:) = y.*weights(i);
end
A = sum(sum(p_fsd));
C = 1/A;
pfsd = p_fsd.*C;

for i = 1:length(alpha)
   label_vec{i} = strcat('$\alpha$',sprintf(' = %g',alpha(i)));
end
label_vec = ["Consolidated","Unconsolidated","Ice edge"];
close all
%%
close all
conFigure(11)
f = figure;
bar3(1:nfsd,pfsd')
yticks('auto')
zlabel('Sea ice area fraction [$\%$]')
xlabel('Ice thickness category')
    hYLabel = get(gca,'XLabel');
    set(hYLabel,'rotation',10,'VerticalAlignment','middle')
ylabel('Floe size category')
    hYLabel = get(gca,'YLabel');
    set(hYLabel,'rotation',-19,'VerticalAlignment','middle')
ylim([0.5,12.5])
view([-140,-180,60])
%legend(label_vec,'Interpreter','latex')
%title('Power law')
exportgraphics(f,'jointfsd.pdf','ContentType','vector')
%%


view([0,-180,0])
    hYLabel = get(gca,'XLabel');
    set(hYLabel,'rotation',0,'VerticalAlignment','middle')
exportgraphics(f,'itd.pdf','ContentType','vector')
    


view([-180,-0,0])
ylabel('Floe size category')
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',0,'VerticalAlignment','middle')
zlabel('Sea ice area fraction [$\%$]')
    exportgraphics(f,'fsd.pdf','ContentType','vector')

% f = figure;
% bar(1:nfsd,sum(pfsd)')
% yticks('auto')
% ylabel('Sea ice fraction')
% xlabel('Floe size category')
%ylim([0.5,12.5])

%%
f = figure;
bar(1:ncat,sum(pfsd'))
yticks('auto')
ylabel('Sea ice fraction')
xlabel('Ice thickness category')
%ylim([0.5,12.5])