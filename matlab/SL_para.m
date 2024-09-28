%% calculate how many particles in the bubble
initialR = 4.5e-6;      % radius m
kB=1.380649e-23;    % boltzman constant
Pinfty = 101325;        % 1 atm
Tinfty = 300;           % K
V=4/3*pi*initialR^3;

%%
N=Pinfty*V/(kB*Tinfty);
Nsim=1e5;
Nratio=N/Nsim
rratio=Nratio^(1/3)