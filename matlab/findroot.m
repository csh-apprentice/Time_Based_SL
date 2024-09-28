function root=findroot(xa,xb,xinit,y,tolerance,countmax)
% find the root of the function xa*x^4+xb*x^2=y
xi=xinit;
count=0;
while abs(xi^4*xa+xi^2*xb-y)>tolerance &  count<countmax
    fx=xi^4*xa+xi^2*xb-y;
    Dfx=4*xi^3*xa+2*xi*xb;
    xi=xi-fx/Dfx;
    count=count+1;
end
root=xi;

end