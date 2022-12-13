function [cp,cl,cdp]=wind_tunnel_fun()
close all
clear
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
parameters = load('setup.mat');
data = load('Data_group_2.mat');
p = data.group_2.p;
Uinf = data.group_2.Uinf;
AoA = data.group_2.AoA;
b = parameters.setup.span;
c = parameters.setup.chord;
rho = parameters.setup.rho;
coord = parameters.setup.coord_taps;
nu = 1.516e-5;
xfoil_cp_0_5 = load('visc_alfa0_velo_5');
xfoil_cp_0_15 = load('visc_alfa0_velo_15');
xfoil_cp_0_25 = load('visc_alfa0_velo_25');
xfoil_cp_5_5 = load('visc_alfa5_velo_5');
xfoil_cp_5_15 = load('visc_alfa5_velo_15');
xfoil_cp_5_25 = load('visc_alfa5_velo_25');
xfoil_cp_10_5 = load('visc_alfa10_velo_5');
xfoil_cp_10_15 = load('visc_alfa10_velo_15');
xfoil_cp_10_25 = load('visc_alfa10_velo_25');
i_xfoil_cp_0 = load('invisc_alfa0');
i_xfoil_cp_5 = load('invisc_alfa5');
i_xfoil_cp_10 = load('invisc_alfa10');
% p=getfield(group_2,'p');
% Uinf=getfield(group_2,'Uinf');
% AoA=getfield(group_2,'AoA');
%%
% Area of airfoil profile 
% Définissez le code Naca 4 digits
code = '0018';

% Déterminez les paramètres du profil d'aile en fonction du code Naca 4 digits
m = str2double(code(1)) / 100; % Épaisseur relative
p_area = str2double(code(2)) / 10; % Position de l'épaisseur relative
t = str2double(code(3:4)) / 100; % Taux de cambrure

% Générez les coordonnées x et y du profil d'aile Naca 4 digits
x = linspace(0, 1, 101);
y = (t/0.2) * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x.^2 + 0.2843 * x.^3 - 0.1015 * x.^4);

% Appliquer la forme de l'épaisseur relative et la position de l'épaisseur relative au profil d'aile
y_c = y .* (1 + 2 .* m .* (p_area - x));

% Ajouter la moitié inférieure du profil d'aile
x = [x, fliplr(x)];
y = [y_c, -fliplr(y_c)];

% Calculez l'aire du profil d'aile en utilisant la formule des trapèzes
n = length(x);
aire = 0;
for i = 1:n-1
    aire = aire + (x(i+1) - x(i)) * (y(i+1) + y(i)) / 2;
end
aire = aire*0.45;
% Affichez l'aire calculée
% disp(aire);
%%
%correct the free stream velocity
K1=0.52;
S=2.5*1.8;%area of the wind tunel section 
VB=aire*b;
eps=K1*VB/S^(3/2);
U=Uinf*(1+eps);
Re = c*U/nu;
pinf=0;%[N/m^2], ambient pressure
cp=zeros(size(p));
%AoAtt = [0 5 10 0 5 10 0 5 10];

alpha=2;
for i=1:size(p,1)
    for j=1:size(p,2)
        cp(i,j)=2*(p(i,j)-pinf)/(rho*U(i)^2);
    end
end

%% Graphes CP
% 0 deg
figure
plot(i_xfoil_cp_0(:,1),i_xfoil_cp_0(:,3),'--');
hold on
plot(coord(1,:)/c,cp(1,:),'o');
plot(xfoil_cp_0_5(:,1),xfoil_cp_0_5(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
title('$C_p$ with $Re=1.2358e5$ for 0$^\circ$AoA');
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
grid on
figure
plot(i_xfoil_cp_0(:,1),i_xfoil_cp_0(:,3),'--');
hold on
plot(coord(1,:)/c,cp(4,:),'o');
plot(xfoil_cp_0_15(:,1),xfoil_cp_0_15(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
title('$C_p$ with $Re=3.7074e5$ for 0$^\circ$AoA');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
figure
plot(i_xfoil_cp_0(:,1),i_xfoil_cp_0(:,3),'--');
hold on
plot(coord(1,:)/c,cp(7,:),'o');
plot(xfoil_cp_0_25(:,1),xfoil_cp_0_25(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
title('$C_p$ with $Re=6.1790e5$ for 0$^\circ$AoA');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
% 5 deg
figure
plot(i_xfoil_cp_5(:,1),i_xfoil_cp_5(:,3),'--');
hold on
plot(coord(1,:)/c,cp(alpha,:),'o');
plot(xfoil_cp_5_5(:,1),xfoil_cp_5_5(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
title('$C_p$ with $Re=1.2358e5$ for 5$^\circ$AoA');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
figure
plot(i_xfoil_cp_5(:,1),i_xfoil_cp_5(:,3),'--');
hold on
plot(coord(1,:)/c,cp(5,:),'o');
plot(xfoil_cp_5_15(:,1),xfoil_cp_5_15(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
title('$C_p$ with $Re=3.7074e5$ for 5$^\circ$AoA');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
figure
plot(i_xfoil_cp_5(:,1),i_xfoil_cp_5(:,3),'--');
hold on
plot(coord(1,:)/c,cp(8,:),'o');
plot(xfoil_cp_5_25(:,1),xfoil_cp_5_25(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
title('$C_p$ with $Re=6.1790e5$ for 5$^\circ$AoA');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
% 10 deg
figure
plot(i_xfoil_cp_10(:,1),i_xfoil_cp_10(:,3),'--');
hold on
plot(coord(1,:)/c,cp(3,:),'o');
plot(xfoil_cp_10_5(:,1),xfoil_cp_10_5(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
title('$C_p$ with $Re=1.2358e5$ for 10$^\circ$AoA');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
figure
plot(i_xfoil_cp_10(:,1),i_xfoil_cp_10(:,3),'--');
hold on
plot(coord(1,:)/c,cp(6,:),'o');
plot(xfoil_cp_10_15(:,1),xfoil_cp_10_15(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
title('$C_p$ with $Re=3.7074e5$ for 10$^\circ$AoA');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
figure
plot(i_xfoil_cp_10(:,1),i_xfoil_cp_10(:,3),'--');
hold on
plot(coord(1,:)/c,cp(9,:),'o');
plot(xfoil_cp_10_25(:,1),xfoil_cp_10_25(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
title('$C_p$ with $Re=6.1790e5$ for 10$^\circ$AoA');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
%% Calcul Cl and Cdp
setup = [0 U(1);5 U(2);10 U(3);...
    0 U(4);5 U(5);10 U(6);...
    0 U(7);5 U(8);10 U(9)];
for j = 1:9
    deltax = zeros(1,size(coord,2));
    deltay = zeros(1,size(coord,2));
    alpha = zeros(1,size(coord,2));
    deltax(1) = (coord(1,2)-c)/2;
    deltax(end) = (c-coord(1,end-1))/2;
    deltay(1) = coord(2,2)/2;
    deltay(end) = -coord(2,end-1);
    alpha(1) = atan2(deltay(1),deltax(1));
    alpha(end) = atan(deltay(end)/deltax(end));
    for i = 2:1:size(coord,2)-1
        deltay(i) = (coord(2,i+1)-coord(2,i-1))/2;
        deltax(i) = (coord(1,i+1)-coord(1,i-1))/2;
        if deltax(i)>=0 && deltay(i)>=0
            alpha(i) = atan(deltay(i)/deltax(i));
        elseif deltax(i)<=0 && deltay(i)>=0
            alpha(i) = atan2(deltay(i),deltax(i));
        elseif deltax(i)<=0 && deltay(i)<=0
            alpha(i) = 2*pi + atan2(deltay(i),deltax(i));
        elseif deltax(i)>=0 && deltay(i)<=0
            alpha(i) = 2*pi + atan(deltay(i)/deltax(i));
        end
    end
    C_l(j) = 2 * (sum(-p(j,:) .* deltax .* cos(alpha -...
        setup(j,1)*pi/180))) / (c*rho*setup(j,2)^2);
    C_dp(j) = 2 * (sum(p(j,:) .* deltay .* sin(alpha -...
        setup(j,1)*pi/180))) / (c*rho*setup(j,2)^2);
end
C_l_xfoil = [0 0.5943 0.9804 0 0.5188 1.0670 0 0.5285 1.0755].';
%% Graphes Cl and CDp
figure
plot(setup(1:3,1),C_l(1:3),'o')
hold on
plot(setup(1:3,1),C_l_xfoil(1:3),'--o')
legend({'Experimental data','Xfoil (viscid)'},'location','northwest');
xlabel('AoA [$^\circ$]')
ylabel('$C_l$ [-]')
title('$Re = 1.2e5$')
grid on
figure
plot(setup(4:6,1),C_l(4:6),'o')
hold on
plot(setup(1:3,1),C_l_xfoil(4:6),'--o')
legend({'Experimental data','Xfoil (viscid)'},'location','northwest');
xlabel('AoA [$^\circ$]')
ylabel('$C_l$ [-]')
title('$Re = 3.7e5$')
grid on
figure
plot(setup(1:3,1),C_l(7:9),'o')
hold on
plot(setup(1:3,1),C_l_xfoil(7:9),'--o')
legend({'Experimental data','Xfoil (viscid)'},'location','northwest');
xlabel('AoA [$^\circ$]')
ylabel('$C_l$ [-]')
title('$Re = 6.0e5$')
grid on

