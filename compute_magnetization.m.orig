function Sz = compute_magnetization(psi_z)
L = size(psi_z,1);
Sz = 0; 
for  i = 1 : L
    for j = 1 : L
        Sz = Sz +psi_z(i,j);
    end
end
Sz = Sz/(L^2);