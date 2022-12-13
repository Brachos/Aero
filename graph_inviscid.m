function []=graph_inviscid(N)
feature('DefaultCharacterSet','UTF8');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex'); 
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultTextFontsize',13);
set(groot, 'defaultAxesFontsize',13);
set(groot, 'defaultLegendFontsize',13);
set(groot, 'defaultLegendLocation','best');
set(0, 'DefaultLineLineWidth', 1.3);
print=0;%saves the figure
pt = 12;
%evolution of the cl wrt the angle of attack and cp wrt x/c.
NACA=[0 0 1 8];
c=1;
V_inf=50;
alpha=[-15 -10 -5 0 5 10 15]*pi/180;
cl_panel=zeros(7,1);
figure
for i=1:7
    [cl,cp,X,N]=my_panel_method(NACA,c,alpha(i),V_inf,N);
    cl_panel(i)=cl;
    if i>=4
        plot(X(1:N)/c,cp);
        hold on
       
    end
end
set(gca, 'YDir','reverse');
grid on
xlabel('x/c [-]');
ylabel('$c_p$ [-]');
legend('$\alpha=0^\circ$','$\alpha=5^\circ$','$\alpha=10^\circ$','$\alpha=15^\circ$');
if print==1
    hgexport(gcf,'cp_panel');
end
hold off
cl_conf=zeros(size(cl_panel));
cl_thin=zeros(size(cl_panel));
for i=1:length(cl_panel)
    cl_conf(i)=2*pi*sin(alpha(i));
    cl_thin(i)=2*pi*sin(alpha(i));
end
figure
plot(alpha*180/pi,cl_panel);
hold on 
plot(alpha*180/pi,cl_conf);
hold on
plot(alpha*180/pi,cl_thin);
xlabel('$\alpha[^\circ]$');
ylabel('$c_l$ [-]');
legend('Panel Method','Conformal Mapping','Thin Airfoil Theory','location','northwest');
grid on
if print==1
    hgexport(gcf,'cl_inviscid');
end
end