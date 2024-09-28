clear;

% Define symbolic variable t
syms t;

% Define normailize paramters
initialR =4.5e-6;      % radius m
Pinfty = 101325;        % 1 atm
Tinfty = 293;           % K
rhoinfty = 1e3;         % water density kg/m^3

R0=initialR;
P0=Pinfty;
T0=Tinfty;
U0=sqrt(P0/rhoinfty);
rho0=rhoinfty;
t0=R0/U0;
k0=Pinfty*U0*R0/Tinfty;
alpha0=U0*R0;
miu0=Pinfty*R0/U0;
sigma0=Pinfty*R0;

% Constants
mygamma = 1.4;          % specific heat ratio
%mygamma= 5/3;           % noble gas

myC = 1481/U0;             % water sound speed m/s
NBC =1.318;            % Hyperparameters here
alpha1 = 0.148e-6/alpha0;      % thermal diffusivity m^2/s
% A = 5.528e-5/(k0/T0);           % air W/mK^2, kg=AT+B
% B = 1.165e-2/k0;           % air W/mK, kg=AT+B 
A = 2.682e-5/(k0/T0);           % argon W/mK^2, kg=AT+B 
B = 1.346e-2/k0;           % argon W/mK, kg=AT+B 
%A = 24.033e-5/(k0/T0);           % Helium W/mK^2, kg=AT+B 
%B = 1.036e-1/k0;           % Helium W/mK, kg=AT+B 

k1 = 0.61/k0;             % water W/mK^2
%k1 = 0.62/k0;             % water W/mK^2
PA = 1.2 * Pinfty/P0;    % ultrasonic amplitude
f = 26.5e3;             % frequency Hz
w = 2 * pi * f;

miu = 0.001/miu0;         % dynamics viscosity Pa*s
%miu = 0.007/miu0;         % dynamics viscosity Pa*s
%miu = 0.000853/miu0;         % dynamics viscosity Pa*s


sigma = 0.072/sigma0;         % surface tension N/m
%sigma = 0.03/sigma0;         % Putterman data surface tension N/m

initialdelta = 0.3 * initialR/R0; % radius m
initialrhog = 1.177/rho0;    % air 300 K 1atm, kg/m^3
%initialrhog =1.603/rho0;     % aragon 300K 1atm, kg/m^3
%initialrhog =0.164/rho0;     % helium 300K 1atm, kg/m^3
Pinitial=Pinfty+2*sigma*sigma0/R0;
a = (4/3) * pi * initialR^3 * initialrhog * 5 / (4 * pi) * (1 - NBC);

shift=0.005;
% Define time parameters
time2run = 10*1e-5/t0;        % simulation time
dt = time2run*1e-6;              % time step

opts = odeset('MaxStep', dt);

yInit = [1.0, 0, Pinitial/P0,1.0,0.3];
interval = [0 time2run];

ySol = ode45(@(t,x) thermo_uniform_ode(t, x, mygamma, myC, Pinfty, Tinfty, rhoinfty, NBC, alpha1, A, B, k1, PA, f, w, miu, sigma, initialR, initialdelta, initialrhog, a,shift,t0), interval, yInit, opts);
NBClist=ySol.y(3,:).*ySol.y(1,:).^3./ySol.y(4,:);

tmax=max(ySol.y(4,:));
figure(1)
clf
plot(ySol.x * 1e6*t0, ySol.y(1,:) *R0* 1e6, 'k')
%plot(ySol.x * 1e6, ySol.y(1,:) , 'k')
title(sprintf('R0 = %1.3f micrometer', 1e6*initialR))
%xlim([16 38])
xlabel('Time (microsecond)')
ylabel('Radius (micrometer)')

f=1;
set(findall(gcf,'-property','FontSize'),'FontSize',9)
set(findall(gcf,'-property','FontName'),'FontName','Arial')

f_sz = [4,2];
set(f, 'PaperUnits', 'inches')
set(f, 'PaperSize', f_sz)
set(f, 'PaperPositionMode', 'manual')
set(f, 'PaperPosition', [0 0 f_sz(1) f_sz(2)])
print(f, '-dpng', sprintf('test.png', i))

%%
figure(2)
clf
plot(ySol.x * 1e6*t0, ySol.y(4,:)*Tinfty, 'k')
title(sprintf('R0 = %1.3f micrometer', 1e6*initialR))
%xlim([4 100])
xlabel('Time (microsecond)')
ylabel('Temperature (K)')

f=1;
set(findall(gcf,'-property','FontSize'),'FontSize',9)
set(findall(gcf,'-property','FontName'),'FontName','Arial')


figure(3)
clf
plot(ySol.x * 1e6*t0, NBClist, 'k')
title(sprintf('R0 = %1.3f micrometer', 1e6*initialR))
%xlim([4 100])
xlabel('Time (microsecond)')
ylabel('NBC')

f=1;
set(findall(gcf,'-property','FontSize'),'FontSize',9)
set(findall(gcf,'-property','FontName'),'FontName','Arial')
saveas(gcf,"./argon/NBC.jpg")

%%
[ddRlist,ddRInterlist]=anlysis_uniform(ySol.x, mygamma, myC, Pinfty, Tinfty, rhoinfty, NBC, alpha1, A, B, k1, PA, f, w, miu, sigma, initialR, initialdelta, initialrhog, a,shift,U0,t0,ySol);
figure(4)
clf
plot(ySol.x * 1e6*t0, ddRlist*U0/t0, 'k')
title(sprintf('R0 = %1.3f micrometer', 1e6*initialR))
%xlim([4 100])
xlabel('Time (microsecond)')
ylabel('ddR')

f=1;
set(findall(gcf,'-property','FontSize'),'FontSize',9)
set(findall(gcf,'-property','FontName'),'FontName','Arial')

Rlist=ySol.y(1,:);
dRlist=ySol.y(2,:);
constlist=ddRlist.*Rlist./dRlist.^2;
