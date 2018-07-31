clear all
x1 = 0.0;
x2 = 2.0;
x3 = -0.175;
x4 = 1.825;

y1 = 0.0;
y2 = 0.175;
y3 = 2.0;
y4 = 2.175;

plot(x1,y1,'*r')
hold on
plot(x2,y2,'*k')
plot(x3,y3,'*b')
plot(x4,y4,'*c')

costheta = x2-x1;
sintheta = y2-y1;

tantheta = sintheta/costheta;

theta = atan2(sintheta,costheta)*180/pi;

length =sqrt(costheta^2 + sintheta^2);

area = length * 2.0;

detJ = area/4.0;

syms zeta1 zeta2

N1 = (1-zeta1)*(1-zeta2)/4;
N2 = (1+zeta1)*(1-zeta2)/4;
N3 = (1-zeta1)*(1+zeta2)/4;
N4 = (1+zeta1)*(1+zeta2)/4;

x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4;

dxdz1 = diff(x,zeta1)
dxdz2 = diff(x,zeta2)
dydz1 = diff(y,zeta1)
dydz2 = diff(y,zeta2)


J = [dxdz1, dxdz2; dydz1, dydz2];
olddetJ = subs(det(J),[zeta1,zeta2],[0.0,0.0])


