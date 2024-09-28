% Do calculus for in Matlab
% Do calculus for in Matlab
%
% Keller Miksis Equation
% As described in https://ieeexplore-ieee-org.clsproxy.library.caltech.edu/document/8091852
% 
clear;
syms Rb(t) Tb0(t) delta(t) rho0 Ps rhor Pb0
syms mygamma myC Pinfty Tinfty rhoinfty NBC alpha1 A B k1 PA f w miu sigma initialR initialdelta initialrhog a shift deltav2 kgv2 Cterm Tadd

dRtmp = diff(Rb,t);
ddRtmp = diff(Rb,t,t);
dddRtmp=diff(Rb,t,t,t);
dTb0tmp=diff(Tb0,t);

m=(4/3) * pi * initialR^3 * initialrhog;
a = m * 5 / (4 * pi) * (1 - NBC);
rhor = a/Rb^3;
c=(m-(4*pi/5)*a)*(3/(4*pi));
rho0 = c/Rb^3;
Ps = -PA*sin(t*w);
 
% Cterm=(1/(20*(mygamma-1)))*((3*mygamma-2)*dRtmp*ddRtmp*Rb+dddRtmp*Rb^2)*((deltav2/k1)*(rho0+5/14*rhor)+(Rb/(2*kgv2))*(rho0+5/21*rhor));
% Tadd=-1/(40*(mygamma-1)*kgv2)*(rho0+5/21*rhor)*((3*mygamma-2)*dRtmp*ddRtmp*Rb^2+dddRtmp*Rb^3)+Cterm;
eta = (Rb/delta)*(k1/B);
Tbll = -B/A*(1 + eta) + B/A*sqrt((1 + eta)^2 + 2*A/B*(Tb0 + A/(2*B)*Tb0^2 + eta*Tinfty));

T0=Tb0+Cterm;
Tll=Tbll+Tadd;

Pb0 = (Pinfty*NBC*initialR^3/Tinfty)*T0/Rb^3; 
Pb1 = Pb0 - 1/2*(rho0 + 1/2*rhor)*ddRtmp*Rb;
PB = Pb1 - 2*sigma/Rb - 4*miu*dRtmp/Rb; 


% pde Tb0
pdeTb0=dTb0tmp == -(3*(mygamma - 1)*Tb0)/Rb*dRtmp - 6*(mygamma - 1)*k1*(Tbll - Tinfty)/(delta*Rb*Pb0);

% pde Ub
Pw=PB +PA*sin((t+Rb/myC)*w) - Pinfty;
dPw=diff(Pw,t);
pdeUb = (1 - dRtmp/myC)*Rb*ddRtmp+ 3/2*dRtmp*dRtmp*(1 - dRtmp/(3*myC)) == ...
        1/rhoinfty*((1 + dRtmp/myC)*Pw +Rb/myC*dPw);

% pde delta
pdedelta=diff(delta,t)==0;
pdedelta=(1+delta/Rb+3/10*(delta/Rb)^2)*diff(delta,t)==6*alpha1/delta-(2*delta/Rb+1/2*(delta/Rb)^2)*dRtmp...
    -delta*(1+delta/(2*Rb)+1/10*(delta/Rb)^2)*1/(Tbll-Tinfty+shift)*diff(Tbll,t);


% Expression
syms dR ddR 


pdeUb = subs(pdeUb, diff(Rb,t,t), ddR);
pdeUb = subs(pdeUb, diff(Rb,t), dR);
pdeUbsolved = solve(pdeUb, ddR);

disp(pdeUbsolved)

syms dTb0
pdeTb0 = subs(pdeTb0, diff(Rb,t,t), ddR);
pdeTb0 = subs(pdeTb0, diff(Rb,t), dR);
pdeTb0 = subs(pdeTb0, diff(Tb0,t), dTb0);
pdeTb0solved=solve(pdeTb0,dTb0);
disp(pdeTb0solved)

syms ddelta
pdedelta = subs(pdedelta, diff(Rb,t,t), ddR);
pdedelta = subs(pdedelta, diff(Rb,t), dR);
pdedelta = subs(pdedelta, diff(Tb0,t), dTb0);
pdedelta = subs(pdedelta, diff(delta,t), ddelta);
pdedeltasolved=solve(pdedelta,ddelta);
disp(pdedeltasolved)