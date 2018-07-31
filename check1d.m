syms x

% B = [1.0, zeta, 0.5*(3.0*zeta^2 - 1.0)];
% 
% dBdz = diff(B,zeta);

x1 = 72;
x2 = 78;

%Legendre
zeta = (2*x - x1 -x2)/(x2-x1);

%Lagrange
% xc = 0.5*(x1+x2);
% zeta = 2.0*(x-xc)/(x2-x1);


%%%%%%%%%%%%%%%%
B = [1.0, zeta, 0.5*(3.0*zeta^2 - 1.0)];

dBdx = diff(B,x);

xb_in = 75;
sigmax = 25;
term1 = 0.5*((x-xb_in)/sigmax)^2;
phi = 1.0*exp(-term1);

for i=1:3
    f = phi*dBdx(1,i);

    DI = int(f, x, x1, x2);
    eval(DI)
end

for i=1:3
   f = B(i)^2;
   mass = int(f, x, x1, x2);
   mass
end

eval(subs(phi,72+5.3238))
