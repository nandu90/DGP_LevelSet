xlen = 1.0;
ylen = 0.5;

xelem = 80;
yelem = 10;

xangle = 5 * pi/180;
yangle = 95 * pi/180;

edge1 = (xlen/xelem)/cos(xangle)
edge2 = (ylen/yelem)/cos(pi/2 - yangle)

xprojection = xlen/xelem;
yprojection = ylen/yelem;

x = zeros(xelem+1);
y = zeros(yelem+1);

fileID = fopen('mesh.dat','w');
%  fprintf(fileID,'%d %d\n',xelem, yelem);

for j=1:yelem+1
    for i=1:xelem+1
       if(i == 1 & j ~= 1)
           y(i,j) = y(i,j-1) + yprojection;
            len = yprojection/cos(pi/2 - yangle);    
            x(i,j) = x(i,j-1) + len*sin(pi/2 - yangle);
       end
       if(i ~= 1)
           x(i,j) = x(i-1,j)+xprojection;
           len = xprojection/cos(xangle);
           y(i,j) = y(i-1,j)+len*sin(xangle);
       end
       fprintf(fileID,'%12.8f %12.8f\n',x(i,j),y(i,j));
       plot(x(i,j),y(i,j),'*')
       hold on
   end
   
end

fclose(fileID);
