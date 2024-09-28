% Author: David Reza Mittelstein (drmittelstein@gmail.com)
% Medical Engineering, California Institute of Technology, 2020

% Keller-Miksis simulations of bubble cavitaiton in response to ultrasound stimulation (batch)

function Ctermlist = anlysis(tlist,mygamma, myC, Pinfty, Tinfty, rhoinfty, NBC, alpha1, A, B, k1, PA, f, w, miu, sigma, initialR, initialdelta, initialrhog, a,c, shift,deltav2,kgv2,ySol)
    % KM_ode(t, x, c, rho_L, p_vap, p_partial, R0, gamma, f, Pus, Patm)
    % t: Time (s)
    % x: [Radius of bubble (m), d/dt Radius of bubble (m/s)]
    % c: Speed of sound in water (m/s)
    % rho_L: Density of water (kg/m3)
    % p_vap: Vapor pressure of liquid (Pa)
    % p_partial: Partial pressure of gas (Pa)
    % R0: Initial bubble radius (m)
    % gamma
    % f: Frequency of US stimulation (Hz)
    % Pus: Pressure (amplitude) of US signal (Pa)
    % Patm: Pressure of atmosphere (Pa)
    % surf_ten: Surface Tension (Pa/m)
    % nu: Viscocity (Pa-s)
    tsize=size(ySol.y,2);
    Ctermlist=zeros(1,tsize);
    for i=1:tsize
        t=tlist(i);
        R = ySol.y(1,i);  % Radius of bubble (m)
        dR = ySol.y(2,i); % d/dt Radius of bubble (m/s)
        ddR=ySol.y(3,i);
        Tb0=ySol.y(4,i);
        delta=ySol.y(5,i);

        dTb0=(Tinfty*k1*R^2*(6*mygamma - 6)*(Tinfty + (B*((k1*R)/(B*delta) + 1))/A - (B*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (Tinfty*k1*R)/(B*delta)))/B)^(1/2))/A))/(NBC*Pinfty*initialR^3*delta) - (dR*Tb0*(3*mygamma - 3))/R;
        rhor = a/R^3;
        rho0 = c/R^3;
        dddR=(myC*rhoinfty*((3*dR^2*(dR/(3*myC) - 1))/2 - ((dR/myC + 1)*(Pinfty - PA*sin(w*(t + R/myC)) + (2*sigma)/R + (4*dR*miu)/R + ddR*R*(a/(4*R^3) + c/(2*R^3)) - (NBC*Pinfty*initialR^3*Tb0)/(Tinfty*R^3)) - (R*((2*dR*sigma)/R^2 - (4*ddR*miu)/R - dR*ddR*(a/(4*R^3) + c/(2*R^3)) + ddR*R*((3*a*dR)/(4*R^4) + (3*c*dR)/(2*R^4)) + (4*dR^2*miu)/R^2 + PA*w*cos(w*(t + R/myC))*(dR/myC + 1) + (NBC*Pinfty*initialR^3*dTb0)/(Tinfty*R^3) - (3*NBC*Pinfty*dR*initialR^3*Tb0)/(Tinfty*R^4)))/myC)/rhoinfty + ddR*R*(dR/myC - 1)))/(R^2*(a/(4*R^3) + c/(2*R^3)));
        Cterm=(1/(20*(mygamma-1)))*((3*mygamma-2)*dR*ddR*R+dddR*R^2)*((deltav2/k1)*(rho0+5/14*rhor)+(R/(2*kgv2))*(rho0+5/21*rhor));
        Ctermlist(1,i)=Cterm;
   
    end

end



% Do calculus for in Matlab
%
% Keller Miksis Equation
% As described in https://ieeexplore-ieee-org.clsproxy.library.caltech.edu/document/8091852
% 
% syms R(t) Pw
% syms c rho sig k eta P0 Pa w R0
% dRtmp = diff(R,t);
% ddRtmp = diff(R,t,t);
% 
% pw = (P0+2*sig/R0)*(R0/R)^(3*k) - 2*sig/R - 4*eta*dRtmp/R - P0 - Pa*sin(w*t);
% dpw = diff(pw, t);
% 
% KM = (...
%     (1-dRtmp/c)*R*ddRtmp + (3/2-dRtmp/(2*c))*dRtmp^2 == ... 
%     (1/rho)*(1+dRtmp/c)*pw + R/(rho*c)*dpw);
% 
% syms dR ddR
% KM = subs(KM, diff(R,t,t), ddR);
% KM = subs(KM, diff(R,t), dR);
% 
% KMsolved = solve(KM, ddR);
% disp(KMsolved)
% 
% Took the output string of above, and pasted it identically into the code
% for out(2) above