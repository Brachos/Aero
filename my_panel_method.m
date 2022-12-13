function [cl,cp,X,N]=my_panel_method(NACA,c,AOA,V_inf,N)
%2D Hess & Smith method to study NACA four-digit airfoils.

%/!\ the AOA has to be given in radians /!\

%get the parameters of the airfoil
tau=(NACA(3)*10+NACA(4))/100;%thickness ratio
eps=NACA(1)/100;             %maximum camber ratio
if NACA(1)==0 && NACA(2)==0  %chordwise position of maximum camber
    p=c/4; %case of a symetric airfoil
else
    p=NACA(2)/10; 
end                          
               

N=1000;                      %number of panels
xi=linspace(0,2*pi,N+1).';
x_coord=c/2*(cos(xi)+1);

%preallocating space for vectors/matrices
T=zeros(N+1,1);
Y_mean=zeros(N+1,1);
dev_Y_mean=zeros(N+1,1);
X=zeros(N+1,1);
Y=zeros(N+1,1);
l_panel=zeros(N,1);
theta_panel=zeros(N,1);
X_mid=zeros(N,1);
Y_mid=zeros(N,1);
A=zeros(N,N);
B=zeros(N,N);
C=zeros(N,N);
D=zeros(N,N);
E=zeros(N,N);
F=zeros(N,N);
G=zeros(N,N);
M_bloc_1=zeros(N,N);
M_bloc_2=zeros(N,1);
M_bloc_3=zeros(1,N);
N_mat=zeros(N+1,1);
Vt=zeros(N,1);
cp=zeros(N,1);
%%
%define the thickness and the camber of the airfoil
for i=1:N+1
    T(i)=10*tau*c*(0.2969*sqrt(x_coord(i)/c)-0.126*x_coord(i)/c-0.3537*(x_coord(i)/c)^2+0.2843*(x_coord(i)/c)^3-0.1015*(x_coord(i)/c)^4);
    if(x_coord(i)/c>=0 && x_coord(i)/c<=p)
        Y_mean(i)=eps*x_coord(i)/p^2*(2*p-x_coord(i)/c);
    elseif(x_coord(i)/c>=p && x_coord(i)/c<=1)
        Y_mean(i)=eps*(c-x_coord(i))/(1-p)^2*(1+x_coord(i)/c-2*p);
    end
end
%%
%coordinates of the upper and lower surfaces
for i=1:N+1
    if(x_coord(i)/c>=0 && x_coord(i)/c<=p)
        dev_Y_mean(i)=2*eps/p*(1-x_coord(i)/(p*c));
    elseif(x_coord(i)/c>=p && x_coord(i)/c<=1)
        dev_Y_mean(i)=-eps/((1-p)^2)*(1+x_coord(i)/c-2*p)+eps*(c-x_coord(i))/((1-p)^2)*1/c;
    end
end
for i=1:N+1
    if i<floor(N/2)+1%lower surface
        X(i)=x_coord(i)+T(i)/2*sin(atan(dev_Y_mean(i)));
        Y(i)=Y_mean(i)-T(i)/2*cos(atan(dev_Y_mean(i)));
    else%upper surface
        X(i)=x_coord(i)-T(i)/2*sin(atan(dev_Y_mean(i)));
        Y(i)=Y_mean(i)+T(i)/2*cos(atan(dev_Y_mean(i)));
    end
end
%%
%characteristics of each panel
for i=1:N
        l_panel(i)=sqrt((X(i+1)-X(i))^2+(Y(i+1)-Y(i))^2); %length of the ith panel
        theta_panel(i)=atan2((Y(i+1)-Y(i)),(X(i+1)-X(i))); %orientation of the ith panel
        X_mid(i)=X(i)+(X(i+1)-X(i))/2;                    %midpoints of the ith panel
        Y_mid(i)=Y(i)+(Y(i+1)-Y(i))/2;
end
%%
%define matrices
for i=1:N
    for j=1:N
        A(i,j)=-(X_mid(i)-X(j))*cos(theta_panel(j))-(Y_mid(i)-Y(j))*sin(theta_panel(j));
        B(i,j)=(X_mid(i)-X(j))^2+(Y_mid(i)-Y(j))^2;
        C(i,j)=sin(theta_panel(i)-theta_panel(j));
        D(i,j)=cos(theta_panel(i)-theta_panel(j));
        E(i,j)=(X_mid(i)-X(j))*sin(theta_panel(j))-(Y_mid(i)-Y(j))*cos(theta_panel(j));
        F(i,j)=log(1+(l_panel(j)^2+2*A(i,j)*l_panel(j))/B(i,j));
        G(i,j)=atan(E(i,j)*l_panel(j)/(A(i,j)*l_panel(j)+B(i,j)));
    end
end


%construction of the matrix M
for i=1:N
    sum_bloc_2=0;
    for j=1:N
        if(i==j)
            M_bloc_1(i,j)=pi;
            sum_bloc_2=sum_bloc_2+0;
        else
            M_bloc_1(i,j)=C(i,j)*F(i,j)/2-D(i,j)*G(i,j);
            sum_bloc_2=sum_bloc_2+D(i,j)*F(i,j)/2+C(i,j)*G(i,j);
        end
    end
    M_bloc_2(i)=sum_bloc_2;
end
sum_bloc_4=0;
for j=1:N   
    M_bloc_3(j)=D(N,j)*F(N,j)/2+C(N,j)*G(N,j)+F(1,j)*D(1,j)/2+C(1,j)*G(1,j);
    sum_bloc_4=sum_bloc_4+(C(N,j)*F(N,j)/2-D(N,j)*G(N,j)+C(1,j)*F(1,j)/2-D(1,j)*G(1,j));

end
M_bloc_4=-sum_bloc_4;
%assembly
M=1/(2*pi)*[M_bloc_1 M_bloc_2;M_bloc_3 M_bloc_4];

%construction of the vector N
for i=1:N+1
    if(i==N+1)
        N_mat(i)=V_inf*(cos(theta_panel(1)-AOA)+cos(theta_panel(N)-AOA));
    else
        N_mat(i)=V_inf*sin(theta_panel(i)-AOA);
    end
end
%%
%vector containing the unknowns: qi and gamma
u=M^(-1)*N_mat;
q=u(1:N);
gamma=u(N+1);

cl=2*gamma/(V_inf*c)* sum(l_panel);

%compute tangential velocity to the panels
for i=1:N
    som=0;
    for j=1:N
        if(j==i)
            som=som+gamma/2;
        else
            som=som-q(j)/(2*pi)*(D(i,j)*F(i,j)/2+C(i,j)*G(i,j))+gamma/(2*pi)*(C(i,j)*F(i,j)/2-D(i,j)*G(i,j));
        end
    end
    Vt(i)=V_inf*cos(theta_panel(i)-AOA)+som;
    cp(i)=1-(Vt(i)./V_inf)^2;
end
%%
% figure
% plot(X,Y);
% hold on
% scatter(X_mid,Y_mid,'.')
% axis([0 1 -0.5 0.5])


