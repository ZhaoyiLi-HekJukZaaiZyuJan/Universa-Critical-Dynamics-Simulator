function Energy = compute_energy(psi_x,psi_y,psi_z,J_z)

L = size(psi_x,1);
Energy = 0;
for  i = 1 : L
    for j = 1 : L
        %sum right:
        i_right = mod2(i+1,L);
        j_right = j;
        Energy = Energy - (psi_x(i,j)*psi_x(i_right,j_right)+psi_y(i,j)*psi_y(i_right,j_right)+psi_z(i,j)*psi_z(i_right,j_right));
            
        %sum left:
        i_left = mod2(i-1,L);
        j_left = j;
        Energy = Energy - (psi_x(i,j)*psi_x(i_left,j_left)+psi_y(i,j)*psi_y(i_left,j_left)+psi_z(i,j)*psi_z(i_left,j_left));

            
        %sum up:
        i_up = i;
        j_up = mod2(j-1,L);
        Energy = Energy - (psi_x(i,j)*psi_x(i_up,j_up)+psi_y(i,j)*psi_y(i_up,j_up)+psi_z(i,j)*psi_z(i_up,j_up));

            
        %sum down:
        i_down = i;
        j_down = mod2(j+1,L);
        Energy = Energy - (psi_x(i,j)*psi_x(i_down,j_down)+psi_y(i,j)*psi_y(i_down,j_down)+psi_z(i,j)*psi_z(i_down,j_down));
            
        %anisotropy:
        Energy = Energy - J_z*(psi_z(i,j)*psi_z(i,j));
    end
end
Energy = Energy/(2*L^2);
