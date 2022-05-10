function [Dx,Dy,Dz] = Derivative(psi_x,psi_y,psi_z,J_z)

    L = size(psi_x,1);
    Dx = zeros(L,L);
    Dy = zeros(L,L);
    Dz = zeros(L,L);
    for i = 1 : L
        for j = 1 : L
            %sum left rigDt:
            
            k = mod2(i+1,L);
            
            Dx(i,j) = Dx(i,j) - ((1+J_z)*psi_y(i,j)*psi_z(k,j)-psi_z(i,j)*psi_y(k,j));
            Dy(i,j) = Dy(i,j) - (psi_z(i,j)*psi_x(k,j)-(1+J_z)*psi_x(i,j)*psi_z(k,j));
            Dz(i,j) = Dz(i,j) - (        psi_x(i,j)*psi_y(k,j)-psi_y(i,j)*psi_x(k,j));
            
            k = mod2(i-1,L);
            Dx(i,j) = Dx(i,j) - ((1+J_z)*psi_y(i,j)*psi_z(k,j)-psi_z(i,j)*psi_y(k,j));
            Dy(i,j) = Dy(i,j) - (psi_z(i,j)*psi_x(k,j)-(1+J_z)*psi_x(i,j)*psi_z(k,j));
            Dz(i,j) = Dz(i,j) - (        psi_x(i,j)*psi_y(k,j)-psi_y(i,j)*psi_x(k,j));
            
            k = mod2(j+1,L);
                Dx(i,j) = Dx(i,j) - ((1+J_z)*psi_y(i,j)*psi_z(i,k)-psi_z(i,j)*psi_y(i,k));
                Dy(i,j) = Dy(i,j) - (psi_z(i,j)*psi_x(i,k)-(1+J_z)*psi_x(i,j)*psi_z(i,k));
                Dz(i,j) = Dz(i,j) - (        psi_x(i,j)*psi_y(i,k)-psi_y(i,j)*psi_x(i,k));
                
            k= mod2(j-1,L);
            
            Dx(i,j) = Dx(i,j) - ((1+J_z)*psi_y(i,j)*psi_z(i,k)-psi_z(i,j)*psi_y(i,k));
                Dy(i,j) = Dy(i,j) - (psi_z(i,j)*psi_x(i,k)-(1+J_z)*psi_x(i,j)*psi_z(i,k));
                Dz(i,j) = Dz(i,j) - (        psi_x(i,j)*psi_y(i,k)-psi_y(i,j)*psi_x(i,k));
                        
        end
    end