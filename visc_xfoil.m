function [x_u,cp_u,x_l,cp_l]=visc_xfoil(NACA,numpanels,AoA,Re)

%inv_xfoil.mupts the viscid x vs cp of the selected airfoil 4 digit 
%at the given AoA, using Xfoil(viscid flow)
%
%INPUT:
%NACA------[string] string containing the 4 digits of the naca airfoil--[-]
%numpanels-[str
% 
%[x_u,cp_u,x_l,cp_l]=visc_xfoil(NACA,numpanels,AoA,Re)
%
%DESCRIPION:
%gives as outing] number of panels to use-----------------------------[-]
%AoA-------[string] angle of attack-------------------------------------[deg]
%Re--------[string] Reynolds number of the flow-------------------------[-]
%
%OUTPUT:
%x_u-----[numpanel/2x1] x/c coordinate of the upper part of the airfoil---[-]
%cp_u----[numpanel/2x1] cp value of the upper part of the airfoil---------[-]
%x_l-----[numpanel/2x1] x/c coordinate of the lower part of the airfoil---[-]
%cp_l----[numpanel/2x1] cp value of the lower part of the airfoil---------[-]
%
%REFERENCE:
%	


%file in which the results will be saved:

saveCp  = 'Cp_R3A4.txt';   % Pressure coeffucuent
saveDTC = 'UDTC_R3A4.txt'; % Delta star, theta, cf

%check if the file exist and in case delete it:
if exist(saveCp,'file')
    delete(saveCp);
end

%input file that will be used in xfoil:
id_in=fopen('xfoil_input.txt','w+');

%airfoil generation:
fprintf(id_in,['NACA ' NACA '\n']);
fprintf(id_in,'PPAR\n');
fprintf(id_in,['N ' numpanels '\n']);
fprintf(id_in,'\n\n');

%find x vs cp:
fprintf(id_in,'OPER\n');
fprintf(id_in,'VISC\n');
fprintf(id_in,[Re '\n']);
fprintf(id_in,'VPAR\n');
fprintf(id_in,'N 7\n');
fprintf(id_in,'\n');
fprintf(id_in,['ALFA ' AoA '\n']);
fprintf(id_in,['!\n']);
fprintf(id_in,['!\n']);
fprintf(id_in,['CPWR ' saveCp '\n']);
fprintf(id_in,['DUMP ' saveDTC '\n']);

fclose(id_in);

%execute xfoil:
cmd='xfoil.exe < xfoil_input.txt';
[status,result]=system(cmd);
delete('xfoil_input.txt');

id=fopen(saveCp,'r');

data=textscan(id,repmat('%f',[1,3]),'HeaderLines',3,...
                  'CollectOutput',1,...
                  'Delimiter','');

fclose(id);

delete(saveCp);

data=data{1};
x_cp=data(:,1);
y_cp=data(:,2);
cp=data(:,3);

%divide extradoss and intradoss datas:
x_u=x_cp(y_cp>=0);
x_l=x_cp(y_cp<0);

cp_u=cp(y_cp>=0);
cp_l=cp(y_cp<0);
end

