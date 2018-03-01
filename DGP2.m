syms zeta1 zeta2 
syms x1 x2 x3 x4
syms y1 y2 y3 y4

%%%%%%%Shape Function for Mapping quads%%%%%%%%%
N1 = (1-zeta1)*(1-zeta2)/4;
N2 = (1+zeta1)*(1-zeta2)/4;
N3 = (1+zeta1)*(1+zeta2)/4;
N4 = (1-zeta1)*(1+zeta2)/4;


%%%Mapping Function%%%
x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
y = N1*y1 + N2*y2 + N3*y3 + N4*y4;

dxdz1 = diff(x,zeta1)
dxdz2 = diff(x,zeta2)
dydz1 = diff(y,zeta1)
dydz2 = diff(y,zeta2)

J = [dxdz1, dxdz2; dydz1, dydz2];
detJ = det(J);

%%%%Test functions for P1 Legendre Polynomials%%%
B = [1, zeta1, zeta2, zeta1*zeta2, (3*zeta1^2-1)/2,(3*zeta2^2-1)/2, (3*zeta1^2-1)*zeta2/2, (3*zeta2^2-1)*zeta1/2 ((3*zeta1^2-1)/2)*(3*zeta2^2-1)/2];

%%%Mass Matrix%%%
M = sym('M%d%d', [9 9]);
for i=1:9
   for j=1:9
       M(i,j) = B(i)*B(j);
   end
end


ccode(M)

subsM = zeros(9);
z1 = -1;
z2 = -1;
for i=1:9
   for j=1:9
       subsM(i,j) = subs(M(i,j),[zeta1,zeta2],[z1,z2]);
   end
end

subsM

subsB = zeros(1,9);
for i=1:9
    subsB(1,i)=subs(B(i),[zeta1,zeta2],[z1,z2]);
end

subsB