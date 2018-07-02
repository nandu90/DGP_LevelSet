syms zeta1 zeta2
syms x y

syms x1 x2 x3 x4
syms y1 y2 y3 y4

N1 = (1-zeta1)*(1-zeta2)/4;
N2 = (1+zeta1)*(1-zeta2)/4;
N3 = (1-zeta1)*(1+zeta2)/4;
N4 = (1+zeta1)*(1+zeta2)/4;

x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4;