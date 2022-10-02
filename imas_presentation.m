addpath functions
clear all
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