function [c_l,c_dp]=wind_tunnel_fun()
clc
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
set(0, 'DefaultLineLineWidth', 1.4);
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
% D?finissez le code Naca 4 digits
code = '0018';

% D?terminez les param?tres du profil d'aile en fonction du code Naca 4 digits
m = str2double(code(1)) / 100; % ?paisseur relative
p_area = str2double(code(2)) / 10; % Position de l'?paisseur relative
t = str2double(code(3:4)) / 100; % Taux de cambrure

% G?n?rez les coordonn?es x et y du profil d'aile Naca 4 digits
x = linspace(0, 1, 101);
y = (t/0.2) * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x.^2 + 0.2843 * x.^3 - 0.1015 * x.^4);

% Appliquer la forme de l'?paisseur relative et la position de l'?paisseur relative au profil d'aile
y_c = y .* (1 + 2 .* m .* (p_area - x));

% Ajouter la moiti? inf?rieure du profil d'aile
x = [x, fliplr(x)];
y = [y_c, -fliplr(y_c)];

% Calculez l'aire du profil d'aile en utilisant la formule des trap?zes
n = length(x);
aire = 0;
for i = 1:n-1
    aire = aire + (x(i+1) - x(i)) * (y(i+1) + y(i)) / 2;
end
aire = aire*0.45;
% Affichez l'aire calcul?e
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
% title('$Re=1.5e5$ for $\alpha = 0^\circ$');
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
grid on
hgexport(gcf,'Cp_alfa0_Re1.eps')
figure
plot(i_xfoil_cp_0(:,1),i_xfoil_cp_0(:,3),'--');
hold on
plot(coord(1,:)/c,cp(4,:),'o');
plot(xfoil_cp_0_15(:,1),xfoil_cp_0_15(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
% title('$Re=4.5e5$ for $\alpha = 0^\circ$');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
hgexport(gcf,'Cp_alfa0_Re4.eps')
figure
plot(i_xfoil_cp_0(:,1),i_xfoil_cp_0(:,3),'--');
hold on
plot(coord(1,:)/c,cp(7,:),'o');
plot(xfoil_cp_0_25(:,1),xfoil_cp_0_25(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
% title('$Re=7.4e5$ for $\alpha = 0^\circ$');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
hgexport(gcf,'Cp_alfa0_Re7.eps')
% 5 deg
figure
plot(i_xfoil_cp_5(:,1),i_xfoil_cp_5(:,3),'--');
hold on
plot(coord(1,:)/c,cp(alpha,:),'o');
plot(xfoil_cp_5_5(:,1),xfoil_cp_5_5(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
% title('$Re=1.5e5$ for $\alpha = 5^\circ$');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
hgexport(gcf,'Cp_alfa5_Re1.eps')
figure
plot(i_xfoil_cp_5(:,1),i_xfoil_cp_5(:,3),'--');
hold on
plot(coord(1,:)/c,cp(5,:),'o');
plot(xfoil_cp_5_15(:,1),xfoil_cp_5_15(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
% title('$Re=4.5e5$ for $\alpha = 5^\circ$');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
hgexport(gcf,'Cp_alfa5_Re4.eps')
figure
plot(i_xfoil_cp_5(:,1),i_xfoil_cp_5(:,3),'--');
hold on
plot(coord(1,:)/c,cp(8,:),'o');
plot(xfoil_cp_5_25(:,1),xfoil_cp_5_25(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
% title('$Re=7.4e5$ for $\alpha = 5^\circ$');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
hgexport(gcf,'Cp_alfa5_Re7.eps')
% 10 deg
figure
plot(i_xfoil_cp_10(:,1),i_xfoil_cp_10(:,3),'--');
hold on
plot(coord(1,:)/c,cp(3,:),'o');
plot(xfoil_cp_10_5(:,1),xfoil_cp_10_5(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
% title('$Re=1.5e5$ for $\alpha = 10^\circ$');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
hgexport(gcf,'Cp_alfa10_Re1.eps')
figure
plot(i_xfoil_cp_10(:,1),i_xfoil_cp_10(:,3),'--');
hold on
plot(coord(1,:)/c,cp(6,:),'o');
plot(xfoil_cp_10_15(:,1),xfoil_cp_10_15(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
% title('$Re=4.5e5$ for $\alpha = 10^\circ$');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
hgexport(gcf,'Cp_alfa10_Re4.eps')
figure
plot(i_xfoil_cp_10(:,1),i_xfoil_cp_10(:,3),'--');
hold on
plot(coord(1,:)/c,cp(9,:),'o');
plot(xfoil_cp_10_25(:,1),xfoil_cp_10_25(:,3));
set(gca, 'YDir','reverse');
legend({'Xfoil (inviscid)','Experimental data','Xfoil (viscid)'})
% title('$Re=7.4e5$ for $\alpha = 10^\circ$');
grid on
xlabel('$x/c$ [-]')
ylabel('$C_p$ [-]')
hgexport(gcf,'Cp_alfa10_Re7.eps')
%% Alternative method
c_n = zeros(1, 9);
c_a = zeros(1, 9);
for j = 19 : -1  : 1
    delta_x = (coord(1, j)- coord(1, j+1))/c;
    delta_y = (coord(2, j)- coord(2, j+1))/c;
    theta = atan2(delta_y, delta_x);
    c_n(:) = c_n(:) + (cp(:, 40-j) - cp(:, j))*delta_x*cos(theta);
    c_a(:) = c_a(:) + (cp(:, 40-j) + cp(:, j))*delta_y*sin(theta);
end
theta = atan2(coord(1, 20), coord(2, 20)); 

c_n(:) = c_n(:) + cp(:, 20)*coord(1, 20)/c * cos(theta);
c_a(:) = c_a(:) + cp(:, 20)*abs(coord(2, 20))/c * sin(theta);
C_l = c_n.*cos(AoA*pi/180) - c_a.*sin(AoA*pi/180);
C_d = c_n.*sin(AoA*pi/180) + c_a.*cos(AoA*pi/180);

%% Calcul Cl and Cdp
p = [p(:, 20:end) p(:, 1:19)];
coord = [coord(:, 20:end) coord(:, 1:19)];

N = p(:, 1)*coord(1, 1) * (-coord(2, 1)/abs(coord(2, 1)));
A = p(:, 1)*coord(2, 1) * (coord(2, 1)/abs(coord(2, 1)));

for j = 2 : length(p(1, :))
    delta_x = abs(coord(1, j)-coord(1, j-1));
    if delta_x ~= 0 
    delta_y = (coord(2, j)-coord(2, j-1)) * (coord(1, j)-coord(1, j-1))/delta_x;
    else
        delta_y = (coord(2, j)-coord(2, j-1));
    end
    if isnan(delta_y)
        delta_y=0;
    end
    N = N + p(:, j)*delta_x * (-coord(2, j)/abs(coord(2, j)));
    A = A + p(:, j)*delta_y * (coord(2, j)/abs(coord(2, j)));
end

%%
N = N.';
A = A.';
c_l = N.*cos(AoA*pi/180)- A.*sin(AoA*pi/180);
c_dp =  A.*cos(AoA*pi/180) + N.*sin(AoA*pi/180);
c_l = c_l./(0.5*rho*(Uinf.^2)*c);
c_dp = c_dp./(0.5*rho*(Uinf.^2)*c);
c_l_xfoil = [0 0.5665 0.9874 0 0.5222 1.0909 0 0.5332 1.0637].';
c_dp_xfoil = [0.00579 0.00869 0.01593 0.00265 0.00437 0.01086 0.00211...
0.00354 0.00891];
c_dp = C_d;
% c_l = C_l;
%% Graphes Cl and CDp
setup = AoA.';

figure
fig1 = plot(setup(1:3,1),c_l(1:3),'--o','color',[0 0.4470 0.7410]);
colour1 = get(fig1,'color');
hold on
p1 = plot(setup(1:3,1),c_l_xfoil(1:3),'-','color',colour1,'Marker',".",'Markersize',15);
fig2 = plot(setup(4:6,1),c_l(4:6),'--o','color',[0.8500 0.3250 0.0980]);
colour2 = get(fig2,'color');
p2 = plot(setup(1:3,1),c_l_xfoil(4:6),'-','color',colour2,'Marker',".",'Markersize',15);
fig3 = plot(setup(1:3,1),c_l(7:9),'--o','color',[0.9290 0.6940 0.1250]);
colour3 = get(fig3,'color');
p3 = plot(setup(1:3,1),c_l_xfoil(7:9),'-','color',colour3,'Marker',".",'Markersize',15);
xlabel('$\alpha$ [$^\circ$]')
ylabel('$c_l$ [-]')
lgd = legend([p1 p2 p3],{'$Re = 1.5e5$','$Re = 4.5e5$','$Re = 7.4e5$'},'location','northwest');
title(lgd,{'Straight line : $\texttt{Xfoil}$','Dashed line : Experimental data'});
grid on
hgexport(gcf,'Cl_vicid.eps');
% CDp
figure
fig1 = plot(setup(1:3,1),c_dp(1:3),'--o','color',[0 0.4470 0.7410]);
colour1 = get(fig1,'color');
hold on
p1 = plot(setup(1:3,1),c_dp_xfoil(1:3),'-','color',colour1,'Marker',".",'Markersize',15);
fig2 = plot(setup(4:6,1),c_dp(4:6),'--o','color',[0.8500 0.3250 0.0980]);
colour2 = get(fig2,'color');
p2 = plot(setup(1:3,1),c_dp_xfoil(4:6),'-','color',colour2,'Marker',".",'Markersize',15);
fig3 = plot(setup(1:3,1),c_dp(7:9),'--o','color',[0.9290 0.6940 0.1250]);
colour3 = get(fig3,'color');
p3 = plot(setup(1:3,1),c_dp_xfoil(7:9),'-','color',colour3,'Marker',".",'Markersize',15);
xlabel('$\alpha$ [$^\circ$]')
ylabel('$c_{dp}$ [-]')
lgd = legend([p1 p2 p3],{'$Re = 1.5e5$','$Re = 4.5e5$','$Re = 7.4e5$'},'location','northwest');
title(lgd,{'Straight line : $\texttt{Xfoil}$','Dashed line : Experimental data'});
grid on
hgexport(gcf,'Cdp_vicid.eps');

%% Graphs for H, delta_star and theta
fid=fopen('Data_Re1.5_alfa0');
dataRe1_alfa0 = textscan(fid,'%f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

fid=fopen('Data_Re1.5_alfa5');
dataRe1_alfa5 = textscan(fid,'%f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

fid=fopen('Data_Re1.5_alfa10');
dataRe1_alfa10 = textscan(fid,'%f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

fid=fopen('Data_Re4.6_alfa0');
dataRe4_alfa0 = textscan(fid,'%f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

fid=fopen('Data_Re4.6_alfa5');
dataRe4_alfa5 = textscan(fid,'%f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

fid=fopen('Data_Re4.6_alfa10');
dataRe4_alfa10 = textscan(fid,'%f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

fid=fopen('Data_Re7.38_alfa0');
dataRe7_alfa0 = textscan(fid,'%f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

fid=fopen('Data_Re7.38_alfa5');
dataRe7_alfa5 = textscan(fid,'%f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);

fid=fopen('Data_Re7.38_alfa10');
dataRe7_alfa10 = textscan(fid,'%f %f %f %f %f %f %f %f','Headerlines',1);
fclose(fid);
ylab = ["$\delta^\star$ [m]" "$\theta$ [m]" "$C_f$ [-]" "H [-]"];
xlab = "$x/c$ [-]";
leg = {'$Re = 1.5e5$','$Re = 4.6e5$', '$Re = 7.4e5$'};
strexp5 = ["delta5.eps" "theta5.eps" "Cf5.eps" "H5.eps"];
strexp10 = ["delta10.eps" "theta10.eps" "Cf10.eps" "H10.eps"];
y_lim = [0 0.05;0 0.014;-0.01 0.07;1 5];
for i=5:8
%     figure
%     plot(dataRe1_alfa0{2}(1:80),dataRe1_alfa0{i}(1:80));
%     hold on
%     plot(dataRe4_alfa0{2}(1:80),dataRe4_alfa0{i}(1:80));
%     plot(dataRe7_alfa0{2}(1:80),dataRe7_alfa0{i}(1:80));
%     if i==8
%         plot(dataRe1_alfa5{2}(1:80),2.59*ones(1,80),'--')
%     end
%     title('$\alpha = 0^\circ$');
%     xlabel(xlab)
%     ylabel(ylab(i-4))
%     legend(leg);
%     grid on
%     hold off
    
    figure
    plot(dataRe1_alfa5{2}(1:80),dataRe1_alfa5{i}(1:80));
    hold on
    plot(dataRe4_alfa5{2}(1:80),dataRe4_alfa5{i}(1:80));
    plot(dataRe7_alfa5{2}(1:80),dataRe7_alfa5{i}(1:80));
    if i==8
        plot(dataRe1_alfa5{2}(1:80),2.59*ones(1,80),'--')
    end
%     title('$\alpha = 5^\circ$');
    xlabel(xlab)
    ylabel(ylab(i-4))
    legend(leg);
    grid on
    ylim(y_lim(i-4,:));
    hgexport(gcf,strexp5(i-4));
    hold off
    
    
    figure
    plot(dataRe1_alfa10{2}(1:80),dataRe1_alfa10{i}(1:80));
    hold on
    plot(dataRe4_alfa10{2}(1:80),dataRe4_alfa10{i}(1:80));
    plot(dataRe7_alfa10{2}(1:80),dataRe7_alfa10{i}(1:80));
    if i==8
        plot(dataRe1_alfa5{2}(1:80),2.59*ones(1,80),'--')
    end
%     title('$\alpha = 10^\circ$');
    xlabel(xlab)
    ylabel(ylab(i-4))
    legend(leg)
    grid on
    ylim(y_lim(i-4,:));
    hgexport(gcf,strexp10(i-4));
    hold off
end
