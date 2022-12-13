clear
data = load('Data_group_8.mat');
condition = load('setup.mat');

p = data.group_8.p;
Uinf = data.group_8.Uinf; 
AoA = data.group_8.AoA;
coor = condition.setup.coord_taps;
rho = condition.setup.rho;
chord = condition.setup.chord;
span = condition.setup.span;

%% Plot data

Figure1=figure(1); clf; set(Figure1,'defaulttextinterpreter','latex');      
plot(coor(1, :)/chord, p(1:5, :),'linewidth',2.0)
ylabel('p [Pa]')
xlabel('x/c [-]')
set(gca,'fontsize',24,'fontname','Times', 'LineWidth',0.5);
legend('\alpha = 0', '\alpha = 4', '\alpha = 8','\alpha = 12', '\alpha = 16', 'location', 'SouthEast')
set(gca,'XTick',[0:0.25:1]);
set(gca,'XMinorTick','off','YMinorTick','off')
set(gcf, 'paperunits', 'inches');
Lx=8; Ly=6;
set(gcf, 'papersize', [Lx Ly]);
set(gcf, 'PaperPosition', [0.01*Lx 0.01*Ly 1.05*Lx 1.02*Ly]);
%hgexport(Figure1,'V_low_p');

Figure2=figure(2); clf; set(Figure2,'defaulttextinterpreter','latex');      
plot(coor(1, :)/chord, p(6:10, :),'linewidth',2.0)
ylabel('p [Pa]')
xlabel('x/c [-]')
set(gca,'fontsize',24,'fontname','Times', 'LineWidth',0.5);
legend('\alpha = 0', '\alpha = 4', '\alpha = 8','\alpha = 12', '\alpha = 16', 'location', 'SouthEast')
set(gca,'XTick',[0:0.25:1]);
set(gca,'XMinorTick','off','YMinorTick','off')
set(gcf, 'paperunits', 'inches');
Lx=8; Ly=6;
set(gcf, 'papersize', [Lx Ly]);
set(gcf, 'PaperPosition', [0.01*Lx 0.01*Ly 1.05*Lx 1.02*Ly]);
%hgexport(Figure2,'V_high_p');
%% Correcting term
S = 2.5*1.8;
K_1 = 0.52;
mu = 18.5*10^(-6);
V_B = 0;
for i = 1 : length(coor(1,:))-1
    V_B = V_B + abs(coor(1, i) - coor(1, i+1))*0.5*(abs(coor(2, i)) + abs(coor(2, i+1)));
end
V_B = V_B*span;

epsilon = K_1*V_B/S^(3/2);
Uinf = Uinf*(1+epsilon);
Re = rho*chord*Uinf/mu;

%% Integration along chord
N_bis = p(:, 20)*abs(coor(1, 20));
A_bis = -p(:, 20)*coor(2, 20);

for j = 21 : 39
    delta_x = coor(1, j)-coor(1, j-1);
    delta_y = -(coor(2, j)-coor(2, j-1)); %% Get the delta_y of the upper surface
    N_bis = N_bis + (p(:, j) - p(:, 40-j))*delta_x;
    A_bis = A_bis +(p(:, 40-j) + p(:, j))*delta_y;
end

%% Circular integration
p = [p(:, 20:end) p(:, 1:19)];
coor = [coor(:, 20:end) coor(:, 1:19)];

N = p(:, 1)*coor(1, 1) * (-coor(2, 1)/abs(coor(2, 1)));
A = p(:, 1)*coor(2, 1) * (coor(2, 1)/abs(coor(2, 1)));

for j = 2 : length(p(1, :))
    delta_x = abs(coor(1, j)-coor(1, j-1));
    delta_y = (coor(2, j)-coor(2, j-1)) * (coor(1, j)-coor(1, j-1))/delta_x;
    N = N + p(:, j)*delta_x * (-coor(2, j)/abs(coor(2, j)));
    A = A + p(:, j)*delta_y * (coor(2, j)/abs(coor(2, j)));
end

%%
c_l = N*cos(AoA*pi/180) - A*sin(AoA*pi/180);
c_p = N*sin(AoA*pi/180) + A*cos(AoA*pi/180);
c_l = c_l/(0.5*rho*Uinf.^2*chord);
c_p = c_p/(0.5*rho*Uinf.^2*chord);

c_l_bis = N_bis*cos(AoA*pi/180) - A_bis*sin(AoA*pi/180);
c_p_bis = N_bis*sin(AoA*pi/180) + A_bis*cos(AoA*pi/180);
c_l_bis = c_l_bis/(0.5*rho*Uinf.^2*chord);
c_p_bis = c_p_bis/(0.5*rho*Uinf.^2*chord);

%% Plot
Figure3=figure(3); clf; set(Figure3,'defaulttextinterpreter','latex'); 
hold on
plot(AoA(1:5), c_p(1:5),'linewidth',2.0)
plot(AoA(6:10), c_p(6:10),'linewidth',2.0)
plot(AoA(1:5), c_p_bis(1:5),'linewidth',2.0)
plot(AoA(6:10), c_p_bis(6:10),'linewidth',2.0)
xlabel('AoA [$^\circ$]')
ylabel('$c_p$ [-]')
set(gca,'fontsize',24,'fontname','Times', 'LineWidth',0.5);
legend('U_{inf} = 5 m/s', 'U_{inf} = 25 m/s', 'U_{inf} = 5 m/s', 'U_{inf} = 25 m/s')
set(gca,'XMinorTick','off','YMinorTick','off')
set(gcf, 'paperunits', 'inches');
Lx=8; Ly=6;
set(gcf, 'papersize', [Lx Ly]);
set(gcf, 'PaperPosition', [0.01*Lx 0.01*Ly 1.05*Lx 1.02*Ly]);
%hgexport(Figure3,'c_p_exp');

Figure4=figure(4); clf; set(Figure4,'defaulttextinterpreter','latex'); 
hold on
plot(AoA(1:5), c_l(1:5),'linewidth',2.0)
plot(AoA(6:10), c_l(6:10),'linewidth',2.0)
plot(AoA(1:5), c_l_bis(1:5),'linewidth',2.0)
plot(AoA(6:10), c_l_bis(6:10),'linewidth',2.0)
xlabel('AoA [$^\circ$]')
ylabel('$c_l$ [-]')
set(gca,'fontsize',24,'fontname','Times', 'LineWidth',0.5);
legend('U_{inf} = 5 m/s', 'U_{inf} = 25 m/s', 'U_{inf} = 5 m/s', 'U_{inf} = 25 m/s')
set(gca,'XMinorTick','off','YMinorTick','off')
set(gcf, 'paperunits', 'inches');
Lx=8; Ly=6;
set(gcf, 'papersize', [Lx Ly]);
set(gcf, 'PaperPosition', [0.01*Lx 0.01*Ly 1.05*Lx 1.02*Ly]);
%hgexport(Figure4,'c_l_exp');