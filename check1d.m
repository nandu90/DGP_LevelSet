format long
syms x

% B = [1.0, zeta, 0.5*(3.0*zeta^2 - 1.0)];
% 
% dBdz = diff(B,zeta);

x1 = 72;
x2 = 78;

zeta = (2*x - x1 -x2)/(x2-x1);

%%%%%%%%%%%%%%%%
%legendre
% B = [1.0, zeta, 0.5*(3.0*zeta^2 - 1.0)];

%Lagrange
B = [zeta*(1-zeta)/2, (1-zeta)*(1+zeta), zeta*(1+zeta)/2];

B= [0.5*(1-zeta),0.5*(1+zeta),0.0]
dBdx = diff(B,x);

xb_in = 75;
sigmax = 25;
term1 = 0.5*((x-xb_in)/sigmax)^2;
phi = 1.0*exp(-term1);

for i=1:2
    f = phi*dBdx(1,i);

    DI = int(f, x, x1, x2);
    eval(DI)
end

for i=1:2
    for j=1:2
        f = B(i)*B(j);
        mass = int(f, x, x1, x2);
        eval(mass)
    end
    end

eval(subs(phi,x,72))
