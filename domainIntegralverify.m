% I have verified the domain integral with the code
% The section of the code that generates the solution
% closest to one coded here is commented out in the
% Assuming the integration over here is exact, the max
% error in between the integrals in this code and that 
% from my domain integral routine is 0.5%
% The code is tested for DGP1 and the first element i.e, (2,2)
% in the code

clear all
format long
xb = 0.6;
yb = 0.6;
rb = 0.2;

x1 = 0.0;
y1 = 0.0;

x2 = 0.01;
y2 = 0.0;

x3 = 0.0;
y3 = 0.01;

x4 = 0.01;
y4 = 0.01;


%%%%%%%Shape Function for Mapping quads%%%%%%%%%
z1 = -1:0.001:1;
z2 = -1:0.001:1;

[zeta1, zeta2] = meshgrid(z1,z2);
func = zeros(size(z1,2),size(z2,2));

for i=1:size(z1,2)
    for j=1:size(z2,2)
    N1 = (1-zeta1(i))*(1-zeta2(j))/4;
    N2 = (1+zeta1(i))*(1-zeta2(j))/4;
    N3 = (1-zeta1(i))*(1+zeta2(j))/4;
    N4 = (1+zeta1(i))*(1+zeta2(j))/4;
    x = N1*x1 + N2*x2 + N3*x3 + N4*x4;
    y = N1*y1 + N2*y2 + N3*y3 + N4*y4;
 
    phi = sqrt((x - xb)^2.0 + (y-yb)^2.0) - rb;

    F1 = 1.0*phi;
    F2 = 0.0*phi;

    gB1 = [0.0, 1.0, 0.0, zeta2(j)];
    gB2 = [0.0, 0.0, 1.0, zeta1(i)];
    
    func1(i,j) = gB1(1)*F1 + gB2(1)*F2;
    func2(i,j) = gB1(2)*F1 + gB2(2)*F2;
    func3(i,j) = gB1(3)*F1 + gB2(3)*F2;
    func4(i,j) = gB1(4)*F1 + gB2(4)*F2;
    end
end

 I1 = trapz(z2,trapz(z1,func1))
 I2 = trapz(z2,trapz(z1,func2))
 I3 = trapz(z2,trapz(z1,func3))
 I4 = trapz(z2,trapz(z1,func4))

% 
% 
% for i=1:4
%     
% end



