a=-7.69e-17;
c=1.9233e-16;
xb=15/(3*a/c+5);
xa=xb*a/c;   %pdf=xa*x^4+xb*x^2

x=linspace(0,1);
y=x.^4*xa+x.^2*xb;

plot(x,y)