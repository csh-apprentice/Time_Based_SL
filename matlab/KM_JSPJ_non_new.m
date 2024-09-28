% Author: David Reza Mittelstein (drmittelstein@gmail.com)
% Medical Engineering, California Institute of Technology, 2020

% Keller-Miksis simulations of bubble cavitaiton in response to ultrasound stimulation (batch)

function out = KM_JSPJ_non_new(t, x, mygamma, myC, Pinfty, Tinfty, rhoinfty, NBC, alpha1, A, B, k1, PA, f, w, miu, sigma, initialR, initialdelta, initialrhog, a,c,shift)
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
    
    R = x(1);  % Radius of bubble (m)
    dR = x(2); % d/dt Radius of bubble (m/s)
    ddR=x(3);
    Tb0=x(4);
    delta=x(5);

    out = zeros(5,1);
    
    out(1) = dR; % d/dt Radius of bubble (m/s) OUT
    out(2)=ddR;
    %out(4)= (Tinfty*k1*R^2*(6*mygamma - 6)*(Tinfty + (B*((k1*R)/(B*delta) + 1))/A - (B*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (Tinfty*k1*R)/(B*delta)))/B)^(1/2))/A))/(NBC*Pinfty*initialR^3*Tb0*delta) - (dR*Tb0*(3*mygamma - 3))/R;
    %out(3)= (Tinfty*k1*R^2*(6*mygamma - 6)*(Tinfty + (B*((k1*R)/(B*delta) + 1))/A - (B*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (Tinfty*k1*R)/(B*delta)))/B)^(1/2))/A))/(NBC*Pinfty*initialR^3*delta) - (dR*Tb0*(3*mygamma - 3))/R;
    out(4)=(Tinfty*k1*R^2*(6*mygamma - 6)*(Tinfty + (B*((k1*R)/(B*delta) + 1))/A - (B*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (Tinfty*k1*R)/(B*delta)))/B)^(1/2))/A))/(NBC*Pinfty*initialR^3*delta) - (dR*Tb0*(3*mygamma - 3))/R;
    dTb0=out(4);


    %out(2) = -((3*dR^2*(dR/(3*myC) - 1))/2 - ((dR/myC + 1)*(Pinfty + (2*sigma)/R - PA*sin(t*w)*(t + R/myC) + (4*dR*miu)/R - (NBC*Pinfty*initialR^3*Tb0)/(Tinfty*R^3)) - (R*((2*dR*sigma)/R^2 + PA*sin(t*w)*(dR/myC + 1) + (4*dR^2*miu)/R^2 + PA*w*cos(t*w)*(t + R/myC) + (NBC*Pinfty*initialR^3*dTb0)/(Tinfty*R^3) - (3*NBC*Pinfty*dR*initialR^3*Tb0)/(Tinfty*R^4)))/myC)/rhoinfty)/(R*(dR/myC - 1) - (R*(dR/myC + 1)*(a/(4*R^3) + (initialR^3*initialrho0)/(2*R^3)) + (R*(dR*(a/(4*R^3) + (initialR^3*initialrho0)/(2*R^3)) - R*((3*a*dR)/(4*R^4) + (3*dR*initialR^3*initialrho0)/(2*R^4)) + (4*miu)/R))/myC)/rhoinfty); % d2/dt2 Radius of bubble(m/s2)
    %out(2)=-((3*dR^2*(dR/(3*myC) - 1))/2 - ((dR/myC + 1)*(Pinfty - PA*sin(w*(t + R/myC)) + (2*sigma)/R + (4*dR*miu)/R - (NBC*Pinfty*initialR^3*Tb0)/(Tinfty*R^3)) - (R*((2*dR*sigma)/R^2 + (4*dR^2*miu)/R^2 + PA*w*cos(w*(t + R/myC))*(dR/myC + 1) + (NBC*Pinfty*initialR^3*dTb0)/(Tinfty*R^3) - (3*NBC*Pinfty*dR*initialR^3*Tb0)/(Tinfty*R^4)))/myC)/rhoinfty)/(R*(dR/myC - 1) - ((R*(dR*(initialrhog/2 + (5*initialR^3*initialrhog*(NBC - 1))/(12*R^3)) + (4*miu)/R - (5*dR*initialR^3*initialrhog*(NBC - 1))/(4*R^3)))/myC + R*(initialrhog/2 + (5*initialR^3*initialrhog*(NBC - 1))/(12*R^3))*(dR/myC + 1))/rhoinfty);
    %out(2)=-((3*dR^2*(dR/(3*myC) - 1))/2 - ((dR/myC + 1)*(Pinfty - PA*sin(w*(t + R/myC)) + (2*sigma)/R + (4*dR*miu)/R - (NBC*Pinfty*initialR^3*Tb0)/(Tinfty*R^3)) - (R*((2*dR*sigma)/R^2 + (4*dR^2*miu)/R^2 + PA*w*cos(w*(t + R/myC))*(dR/myC + 1) + (NBC*Pinfty*initialR^3*dTb0)/(Tinfty*R^3) - (3*NBC*Pinfty*dR*initialR^3*Tb0)/(Tinfty*R^4)))/myC)/rhoinfty)/(R*(dR/myC - 1) - ((R*(dR*(initialrhog/2 + (5*initialR^3*initialrhog*(NBC - 1))/(12*R^3)) + (4*miu)/R - (5*dR*initialR^3*initialrhog*(NBC - 1))/(4*R^3)))/myC + R*(initialrhog/2 + (5*initialR^3*initialrhog*(NBC - 1))/(12*R^3))*(dR/myC + 1))/rhoinfty);
    out(3)=(myC*rhoinfty*((3*dR^2*(dR/(3*myC) - 1))/2 - ((dR/myC + 1)*(Pinfty - PA*sin(w*(t + R/myC)) + (2*sigma)/R + (4*dR*miu)/R + ddR*R*(a/(4*R^3) + c/(2*R^3)) - (NBC*Pinfty*initialR^3*Tb0)/(Tinfty*R^3)) - (R*((2*dR*sigma)/R^2 - (4*ddR*miu)/R - dR*ddR*(a/(4*R^3) + c/(2*R^3)) + ddR*R*((3*a*dR)/(4*R^4) + (3*c*dR)/(2*R^4)) + (4*dR^2*miu)/R^2 + PA*w*cos(w*(t + R/myC))*(dR/myC + 1) + (NBC*Pinfty*initialR^3*dTb0)/(Tinfty*R^3) - (3*NBC*Pinfty*dR*initialR^3*Tb0)/(Tinfty*R^4)))/myC)/rhoinfty + ddR*R*(dR/myC - 1)))/(R^2*(a/(4*R^3) + c/(2*R^3)));
 
    

    eta = (R/delta)*(k1/B);
    Tbll = -B/A*(1 + eta) + B/A*sqrt((1 + eta)^2 + 2*A/B*(Tb0 + A/(2*B)*Tb0^2 + eta*Tinfty));
    
    %out(5)=0;
    if (abs(Tbll-Tinfty)<=shift)
        out(5)=0;
    else 
        %out(4)=-(dR*((2*delta)/R + delta^2/(2*R^2)) - (6*alpha1)/delta + (delta*((dR*k1)/(A*delta) - (B*((2*A*(dTb0 + (A*dTb0*Tb0)/B + (Tinfty*dR*k1)/(B*delta)))/B + (2*dR*k1*((k1*R)/(B*delta) + 1))/(B*delta)))/(2*A*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (Tinfty*k1*R)/(B*delta)))/B)^(1/2)))*(delta/(2*R) + delta^2/(10*R^2) + 1))/(Tinfty + (B*((k1*R)/(B*delta) + 1))/A - (B*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (Tinfty*k1*R)/(B*delta)))/B)^(1/2))/A))/(delta/R + (3*delta^2)/(10*R^2) + (delta*((B*((2*k1*R*((k1*R)/(B*delta) + 1))/(B*delta^2) + (2*A*Tinfty*k1*R)/(B^2*delta^2)))/(2*A*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (Tinfty*k1*R)/(B*delta)))/B)^(1/2)) - (k1*R)/(A*delta^2))*(delta/(2*R) + delta^2/(10*R^2) + 1))/(Tinfty + (B*((k1*R)/(B*delta) + 1))/A - (B*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (Tinfty*k1*R)/(B*delta)))/B)^(1/2))/A) + 1);
        out(5)=-(dR*((2*delta)/R + delta^2/(2*R^2)) - (6*alpha1)/delta + (delta*((dR*k1)/(A*delta) - (B*((2*A*(dTb0 + (A*dTb0*Tb0)/B + (Tinfty*dR*k1)/(B*delta)))/B + (2*dR*k1*((k1*R)/(B*delta) + 1))/(B*delta)))/(2*A*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (Tinfty*k1*R)/(B*delta)))/B)^(1/2)))*(delta/(2*R) + delta^2/(10*R^2) + 1))/(Tinfty + (B*((k1*R)/(B*delta) + 1))/A - (B*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (Tinfty*k1*R)/(B*delta)))/B)^(1/2))/A))/(delta/R + (3*delta^2)/(10*R^2) + (delta*((B*((2*k1*R*((k1*R)/(B*delta) + 1))/(B*delta^2) + (2*A*Tinfty*k1*R)/(B^2*delta^2)))/(2*A*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (Tinfty*k1*R)/(B*delta)))/B)^(1/2)) - (k1*R)/(A*delta^2))*(delta/(2*R) + delta^2/(10*R^2) + 1))/(Tinfty + (B*((k1*R)/(B*delta) + 1))/A - (B*(((k1*R)/(B*delta) + 1)^2 + (2*A*(Tb0 + (A*Tb0^2)/(2*B) + (Tinfty*k1*R)/(B*delta)))/B)^(1/2))/A) + 1);
        %out(5)=0;
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