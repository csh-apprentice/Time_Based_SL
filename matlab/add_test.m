clear;

%% Parameters when generating atoms
N=1e5; % atoms to generate
a=-7.69e-17;
c=1.9233e-16;
xb=15/(3*a/c+5);
xa=xb*a/c;   %pdf=xa*x^4+xb*x^2
Rb=4.5e-6;
N=256*2; % number of intervals in the Ziggurat algorithm
Natoms=1e5; % number of atoms we sample
xlist=zeros(N+1,1);
ylist=zeros(N+1,1);   % pdf list regarding to the xlist
xlist(1)=0.0;
ylist(1)=0.0
%A=1/N; %Area per interval
tolerance=1/N^2;
countmax=50;
error=1.0;
%% generate the xlist
xhigh=2*findroot(xa,xb,0.05,1/N,tolerance,countmax); %an overestimation on the x1
xlow=0.0; % an underestimation on the x1
xguess=(xhigh+xlow)/2;
count=0;
% bisection searching 

%%
while count<countmax
    xlist(2)=xguess;
    ylist(2)=xguess^4*xa+xguess^2*xb;
    Aguess=(1-xlist(2))*ylist(2)+xguess^5/5*xa+xguess^3/3*xb;
    for i=3:N+1
        ylist(i)=ylist(i-1)+Aguess/(1-xlist(i-1));
        xlist(i)=findroot(xa,xb,xlist(i-1),ylist(i),tolerance,countmax);
        if(ylist(i)>xa+xb)
            xhigh=xguess;
            xguess=(xlow+xhigh)/2;
            ylist(N+1)=xa+xb+1;   % break the outeer loop
            break;
        end
    end   
    error=abs(ylist(N+1)-xa-xb);
    if error<tolerance
        break;
    elseif ylist(N+1)<xa+xb
        xlow=xguess;
        xguess=(xlow+xhigh)/2;
    end
    count=count+1;
end

if count==countmax
    disp("Max Try reached in bisection searching");
end
% 
xlist(2)=xguess;
ylist(2)=xguess^4*xa+xguess^2*xb;
Aguess=(1-xlist(2))*ylist(2)+xguess^5/5*xa+xguess^3/3*xb;
for i=3:N+1
    ylist(i)=ylist(i-1)+Aguess/(1-xlist(i-1));
    xlist(i)=findroot(xa,xb,xlist(i-1),ylist(i),tolerance,countmax);
end

% for i=N:-1:2`
%     Ai=(i-1)*A;
%     xi=xlist(i+1);
%     count=0;
%     while abs(xi^5/5*xa+xi^3/3*xb-Ai)>tolerance &  count<countmax
%         fx=xi^5/5*xa+xi^3/3*xb-Ai;
%         Dfx=xi^4*xa++xi^2*xb;
%         xi=xi-fx/Dfx;
%         count=count+1;
%     end
%     xlist(i)=xi;
%     ylist(i)=xa*xi^4+xb*xi^2;
% end

[rlist,plist]=add_atoms(xlist,ylist,xa,xb,Rb,N,Natoms);
rmax=max(rlist)
h=histogram(rlist,100)
%scatter3(plist(:,1),plist(:,2),plist(:,3));
%%
density1=sum(rlist<0.2)/0.2^3;
density2=sum(rlist<0.3)/0.3^3;
density3=sum(rlist<0.9)/0.9^3;
sum(rlist<0.05)