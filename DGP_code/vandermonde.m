vand = [1, -1, -1, 1;
        1, 1, -1, -1;
        1, -1, 1, -1;
        1, 1, 1, 1;];
    
phig = [0.0; 0.1253; 0.0; 0.1253];

phi = inv(vand)*phig;

%reconstruct now
recphi = vand*phi;
recphi
