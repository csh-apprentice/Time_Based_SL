clear;
syms Rb(t) Tb0(t) delta(t) Pb0(t) rho0 Ps rhor 
syms mygamma myC Pinfty Tinfty rhoinfty NBC alpha1 A B k1 PA f w miu sigma initialR initialdelta initialrhog a shift t0

dRtmp = diff(Rb,t);
ddRtmp = diff(Rb,t,t);
dTb0tmp=diff(Tb0,t);
dPb0tmp=diff(Pb0,t);

Pb1 = Pb0;
PB = Pb1 - 2*sigma/Rb - 4*miu*dRtmp/Rb; 

eta = (Rb/delta)*(k1/B);
Tbll = -B/A*(1 + eta) + B/A*sqrt((1 + eta)^2 + 2*A/B*(Tb0 + A/(2*B)*Tb0^2 + eta));

% pde Tb0
pdeTb0=dTb0tmp == -(3*(mygamma - 1)*Tb0)*dRtmp/Rb - 6*(mygamma - 1)*k1*(Tbll - 1)/(delta*Rb*Pb0);
% pde Ub
Pw=PB +PA*sin((t+Rb/myC)*w*t0) - 1;
dPw=diff(Pw,t);
pdeUb = (1 - dRtmp/myC)*Rb*ddRtmp+ 3/2*dRtmp*dRtmp*(1 - dRtmp/(3*myC)) == ...
        ((1 + dRtmp/myC)*Pw +Rb/myC*dPw);

% pde delta
pdedelta=diff(delta,t)==0;
pdedelta=(1+delta/Rb+3/10*(delta/Rb)^2)*diff(delta,t)==6*alpha1/delta-(2*delta/Rb+1/2*(delta/Rb)^2)*dRtmp...
    -delta*(1+delta/(2*Rb)+1/10*(delta/Rb)^2)*1/(Tbll-1)*diff(Tbll,t);

% pde Pb0
pdePb0=dPb0tmp==(-3*mygamma*Pb0/Rb)*dRtmp-6*(mygamma - 1)*k1*(Tbll - 1)/(delta*Rb);

% Expression
syms dR ddR dddR
syms dPb0
pdeUb= subs(pdeUb,diff(Rb,t,t,t),dddR);
pdeUb = subs(pdeUb, diff(Rb,t,t), ddR);
pdeUb = subs(pdeUb, diff(Rb,t), dR);
pdeUb = subs(pdeUb, diff(Pb0,t), dPb0);
pdeUbsolved = solve(pdeUb, ddR);

disp(pdeUbsolved)

syms dTb0
pdeTb0 = subs(pdeTb0, diff(Rb,t,t), ddR);
pdeTb0 = subs(pdeTb0, diff(Rb,t), dR);
pdeTb0 = subs(pdeTb0, diff(Tb0,t), dTb0);
pdeTb0solved=solve(pdeTb0,dTb0);
disp(pdeTb0solved)


pdePb0=subs(pdePb0,diff(Rb,t,t), ddR);
pdePb0=subs(pdePb0,diff(Rb,t), dR);
pdePb0 = subs(pdePb0, diff(Pb0,t), dPb0);
pdePb0solved=solve(pdePb0,dPb0);
disp(pdePb0solved)

syms ddelta
pdedelta = subs(pdedelta, diff(Rb,t,t), ddR);
pdedelta = subs(pdedelta, diff(Rb,t), dR);
pdedelta = subs(pdedelta, diff(Tb0,t), dTb0);
pdedelta = subs(pdedelta, diff(delta,t), ddelta);
pdedeltasolved=solve(pdedelta,ddelta);
disp(pdedeltasolved)