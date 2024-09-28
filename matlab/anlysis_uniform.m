% Author: David Reza Mittelstein (drmittelstein@gmail.com)
% Medical Engineering, California Institute of Technology, 2020

% Keller-Miksis simulations of bubble cavitaiton in response to ultrasound stimulation (batch)

function [ddRlist,ddRInterlist] = anlysis_uniform(tlist, mygamma, myC, Pinfty, Tinfty, rhoinfty, NBC, alpha1, A, B, k1, PA, f, w, miu, sigma, initialR, initialdelta, initialrhog, a,shift,U0,t0,ySol)
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
    ddRlist=zeros(1,tsize);
    ddRInterlist=zeros(1,tsize);
    a0=U0/t0;
    Rlist= ySol.y(1,:);
    dRlist= ySol.y(2,:);
    for i=1:tsize
        t=tlist(i);
        R = ySol.y(1,i);  % Radius of bubble (m)
        dR = ySol.y(2,i); % d/dt Radius of bubble (m/s)
        Pb0=ySol.y(3,i);
        Tb0=ySol.y(4,i);
        delta=ySol.y(5,i);
        dPb0=(k1*(6*mygamma - 6)*((B*((k1*R)/(B*delta) + 1))/A - (B*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (k1*R)/(B*delta)))/B)^(1/2))/A + 1))/(R*delta) - (3*dR*mygamma*Pb0)/R;
        ddR=-((3*dR^2*(dR/(3*myC) - 1))/2 - (dR/myC + 1)*((2*sigma)/R - PA*sin(t0*w*(t + R/myC)) - Pb0 + (4*dR*miu)/R + 1) + (R*(dPb0 + (2*dR*sigma)/R^2 + (4*dR^2*miu)/R^2 + PA*t0*w*cos(t0*w*(t + R/myC))*(dR/myC + 1)))/myC)/(R*(dR/myC - 1) - (4*miu)/myC);
        ddRlist(1,i)=ddR;
        if(i<tsize && i>1)
            ddRInterlist(1,i)=(dRlist(i+1)-dRlist(i-1))/(tlist(i+1)-tlist(i-1));
        elseif i==1
            ddRInterlist(1,i)=(dRlist(i+1)-dRlist(i))/(tlist(i+1)-tlist(i));
        else
            ddRInterlist(1,i)=(dRlist(i)-dRlist(i-1))/(tlist(i)-tlist(i-1));
        end
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