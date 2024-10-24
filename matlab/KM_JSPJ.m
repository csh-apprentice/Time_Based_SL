clear;

% Define symbolic variable t
syms t;

% Constants
%mygamma = 1.4;          % specific heat ratio
mygamma= 5/3;           % noble gas
myC = 1481;             % water sound speed m/s
Pinfty = 101325;        % 1 atm
Tinfty = 300;           % K
rhoinfty = 1e3;         % water density kg/m^3
NBC =1.316;            % Hyperparameters here
alpha1 = 0.148e-6;      % thermal diffusivity m^2/s
% A = 5.528e-5;           % air W/mK^2, kg=AT+B
% B = 1.165e-2;           % air W/mK, kg=AT+B
A = 2.682e-5       % aragon W/mK^2, kg=AT+B
B = 1.346e-2;           % aragon W/mK, kg=AT+B 
% A = 24.033e-5;           % Helium W/mK^2, kg=AT+B 
% B = 1.036e-1;           % Helium W/mK, kg=AT+B 
k1 = 0.61;             % water W/mK^2
PA = 1.2* Pinfty;    % ultrasonic amplitude
f = 26.5e3;             % frequency Hz
w = 2 * pi * f;
miu = 1e-6;         % dynamics viscosity Pa*s
miu=0.001;
sigma = 0.0720;         % surface tension N/m
initialR = 4.5e-6;      % radius m
initialdelta = 0.3 * initialR; % radius m
% initialrhog = 1.177;    % air 300K 1atm, kg/m^3
initialrhog = 1.603;    % aragon 300 K 1atm, kg/m^3
% initialrhog =0.164;     % helium 300K 1atm, kg/m^3
Pinitial=Pinfty+2*sigma/initialR;
m=(4/3) * pi * initialR^3 * initialrhog;
a = (4/3) * pi * initialR^3 * initialrhog * 5 / (4 * pi) * (1 - NBC);
c=(m-(4*pi/5)*a)*(3/(4*pi));
shift=0.005;
deltav2=0.1e-6;
kgv2=5; % air, argon
% kgv2=10.0; % helium
% Define time parameters
time2run = 1e-4;        % simulation time
dt = 1e-10;              % time step

opts = odeset('MaxStep', dt);

%%
Pb0initial=Pinfty*NBC;
intialrho0=c/initialR^3;
intialrhor=a/initialR^3;
initialddR=-(Pinfty + (2*sigma)/initialR - (NBC*Pinfty*initialR^3*Tinfty)/(Tinfty*initialR^3))/(initialR*(a/(4*initialR^3) + c/(2*initialR^3)))
initialddR=-(Pinfty+2*(sigma/initialR)-Pb0initial)/(1/2*(intialrho0+1/2*intialrhor)*initialR)
%yInit = [initialR, 0,Tinfty,initialdelta];
%%
%yInit = [initialR, 0,Tinfty,initialdelta];
yInit = [initialR, 0,initialddR,Tinfty,initialdelta];
interval = [0 time2run];

%ySol = ode45(@(t,x) KM_JSPJ_non(t, x, mygamma, myC, Pinfty, Tinfty, rhoinfty, NBC, alpha1, A, B, k1, PA, f, w, miu, sigma, initialR, initialdelta, initialrhog, a,c, shift), interval, yInit, opts);
ySol = ode15s(@(t,x) KM_JSPJ_non_new(t, x, mygamma, myC, Pinfty, Tinfty, rhoinfty, NBC, alpha1, A, B, k1, PA, f, w, miu, sigma, initialR, initialdelta, initialrhog, a,c, shift), interval, yInit, opts);
%ySol = ode45(@(t,x) KM_JSPJ_non_add(t, x, mygamma, myC, Pinfty, Tinfty, rhoinfty, NBC, alpha1, A, B, k1, PA, f, w, miu, sigma, initialR, initialdelta, initialrhog, a, shift,deltav2,kgv2), interval, yInit, opts);


tmax=max(ySol.y(4,:));
%%
% figure(1)
% clf
% plot(ySol.x * 1e6, ySol.y(1,:) * 1e6, 'k')
% title(sprintf('R0 = %1.3f micrometer', 1e6*initialR))
% %xlim([4 100])
% xlabel('Time (microsecond)')
% ylabel('Radius (micrometer)')
% 
% f=1;
% set(findall(gcf,'-property','FontSize'),'FontSize',9)
% set(findall(gcf,'-property','FontName'),'FontName','Arial')
% 
% f_sz = [4,2];
% set(f, 'PaperUnits', 'inches')
% set(f, 'PaperSize', f_sz)
% set(f, 'PaperPositionMode', 'manual')
% set(f, 'PaperPosition', [0 0 f_sz(1) f_sz(2)])
% print(f, '-dpng', sprintf('test.png', i))
% % 
% figure(2)
% clf
% plot(ySol.x * 1e6, ySol.y(4,:), 'k')
% title(sprintf('R0 = %1.3f micrometer', 1e6*initialR))
% %xlim([4 100])
% xlabel('Time (microsecond)')
% ylabel('Uniform Temperature (K)')
% 
% f=1;
% set(findall(gcf,'-property','FontSize'),'FontSize',9)
% set(findall(gcf,'-property','FontName'),'FontName','Arial')

% 
% figure(3)
% clf
% plot(ySol.x * 1e6, ySol.y(2,:), 'k')
% title(sprintf('R0 = %1.3f micrometer', 1e6*initialR))
% %xlim([4 100])
% xlabel('Time (microsecond)')
% ylabel('dR (m/s)')
% 
% f=1;
% set(findall(gcf,'-property','FontSize'),'FontSize',9)
% set(findall(gcf,'-property','FontName'),'FontName','Arial')
%%
Ctermlist=anlysis(ySol.x,mygamma, myC, Pinfty, Tinfty, rhoinfty, NBC, alpha1, A, B, k1, PA, f, w, miu, sigma, initialR, initialdelta, initialrhog, a,c, shift,deltav2,kgv2,ySol);
TempSum=Ctermlist+ySol.y(4,:);
%Cterm=(1/(20*(mygamma-1)))*((3*mygamma-2)*dRlist.*ddRlist.*Rlist+dddRlist.*dRlist.^2).*((deltav2/k1)*(rho0list+5/14*rhorlist)+(Rlist/(2*kgv2)).*(rho0list+5/21*rhorlist));
% figure(4)
% clf
% plot(ySol.x * 1e6, ySol.y(3,:), 'k')
% title(sprintf('R0 = %1.3f micrometer', 1e6*initialR))
% %xlim([4 100])
% xlabel('Time (microsecond)')
% ylabel('ddR')
% 
% f=1;
% set(findall(gcf,'-property','FontSize'),'FontSize',9)
% set(findall(gcf,'-property','FontName'),'FontName','Arial')


% figure(4)
% clf
% plot(ySol.x * 1e6, Ctermlist, 'k')
% title(sprintf('R0 = %1.3f micrometer', 1e6*initialR))
% %xlim([4 100])
% xlabel('Time (microsecond)')
% ylabel('Cterm (K)')
% 
% f=1;
% set(findall(gcf,'-property','FontSize'),'FontSize',9)
% set(findall(gcf,'-property','FontName'),'FontName','Arial')
% 
% figure(5)
% clf
% plot(ySol.x * 1e6, TempSum, 'k')
% title(sprintf('R0 = %1.3f micrometer', 1e6*initialR))
% %xlim([4 100])
% xlabel('Time (microsecond)')
% ylabel('TempSum (K)')
% 
% f=1;
% set(findall(gcf,'-property','FontSize'),'FontSize',9)
% set(findall(gcf,'-property','FontName'),'FontName','Arial')
%%
[ddRmax,tddRmax]=max(ySol.y(3,:))
[Tempsummax,tTempsummax]=max(TempSum)
[Rmax,tRmax]=max(ySol.y(1,:))
[Rmin,tRmin]=min(ySol.y(1,:))
[dRmax,tRmax]=min(ySol.y(2,:))


%[ddRmax,tddRmax]=max(ySol.y(3,:))
[TempMD,tTempMD]=min(abs(TempSum-Tinfty))

TempSum(tRmin-50)

%% Collecting MD data
MDstart=tRmin-50;
MDdt = 1/4*1e-15;              % time step, set to 1/4 fs in MD
MDtime2run = 100*dt+MDdt*3;        % simulation time
MDopts = odeset('MaxStep', MDdt);

%yInit = [initialR, 0,Tinfty,initialdelta];
MDyInit = [ySol.y(1,MDstart), ySol.y(2,MDstart),ySol.y(3,MDstart),ySol.y(4,MDstart),ySol.y(5,MDstart)];
MDinterval = [0 MDtime2run];
MDySol = ode15s(@(t,x) KM_JSPJ_non_new(t, x, mygamma, myC, Pinfty, Tinfty, rhoinfty, NBC, alpha1, A, B, k1, PA, f, w, miu, sigma, initialR, initialdelta, initialrhog, a,c, shift), MDinterval, MDyInit, MDopts);

%%
% figure(6)
% clf
% plot(MDySol.x * 1e6, MDySol.y(1,:) * 1e6, 'k')
% title(sprintf('R0 = %1.3f micrometer', 1e6*initialR))
% %xlim([4 100])
% xlabel('Time (microsecond)')
% ylabel('Radius (micrometer)')

% f=1;
% set(findall(gcf,'-property','FontSize'),'FontSize',9)
% set(findall(gcf,'-property','FontName'),'FontName','Arial')

% f_sz = [4,2];
% set(f, 'PaperUnits', 'inches')
% set(f, 'PaperSize', f_sz)
% set(f, 'PaperPositionMode', 'manual')
% set(f, 'PaperPosition', [0 0 f_sz(1) f_sz(2)])
% print(f, '-dpng', sprintf('test.png', i))


% figure(7)
% clf
% plot(MDySol.x * 1e6, MDySol.y(4,:), 'k')
% title(sprintf('R0 = %1.3f micrometer', 1e6*initialR))
% %xlim([4 100])
% xlabel('Time (microsecond)')
% ylabel('Temperature (K)')

% f=1;
% set(findall(gcf,'-property','FontSize'),'FontSize',9)
% set(findall(gcf,'-property','FontName'),'FontName','Arial')
%% Store the data
pathR='/home/shihan/projects/sonoluminescence/lammps/build/KM_ar/R_data/1_4fs.txt';
pathdR='/home/shihan/projects/sonoluminescence/lammps/build/KM_ar/dR_data/1_4fs.txt';
pathTbl='/home/shihan/projects/sonoluminescence/lammps/build/KM_ar/Tbl_data/1_4fs.txt';

etalist=(MDySol.y(1,:)./MDySol.y(5,:))*(k1/B);
Tblllist = -B/A*(1 + etalist) + B/A*sqrt((1 + etalist).*(1 + etalist) + 2*A/B*(MDySol.y(4,:) + A/(2*B)*MDySol.y(4,:).*MDySol.y(4,:) + etalist*Tinfty));

writematrix(MDySol.y(1,:).'/1e-10,pathR);
writematrix(MDySol.y(2,:).'/1e5,pathdR);
writematrix(Tblllist',pathTbl);

%%

index=MDstart;
delta0=ySol.y(5,index)*1e10
Rb0=ySol.y(1,index)*1e10
Tb0=ySol.y(4,index)
eta=(Rb0/delta0)*(k1/B);


Tbll = -B/A*(1 + eta) + B/A*sqrt((1 + eta)^2 + 2*A/B*(Tb0 + A/(2*B)*Tb0^2 + eta*Tinfty))
fracr=1.0;
%-2*eta*(A/B)*(Tbll-Tinfty)*fracr^2


Tbr=B/A*(-1+sqrt( (1+A/B*Tb0)^2-2*eta*A/B*(Tbll-Tinfty)*fracr^2 ))


%%
