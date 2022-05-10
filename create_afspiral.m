function [psi_x,psi_y,psi_z] = create_afspiral(L,qx,qy,theta,S)

psi_x = zeros(L,L);
psi_y = zeros(L,L);
psi_z = zeros(L,L);

for i = 1 : L
    for j = 1 : L
        Omega_x = normrnd(0,1/sqrt(2*S));
        Omega_y = normrnd(0,1/sqrt(2*S));
        Omega_z = 1;
        
        if mod(i+j,2)
            theta2 = theta;
        else
            theta2 = theta+pi();
        end

        phi = qx*i+qy*j;
        Rot = [cos(theta2)*cos(phi) -sin(phi) sin(theta2)*cos(phi);cos(theta2)*sin(phi) cos(phi) sin(theta2)*sin(phi);-sin(phi) 0 cos(theta2)];
        psi_x(i,j)= Rot(1,:)*[Omega_x ; Omega_y ; Omega_z];  
        psi_y(i,j)= Rot(2,:)*[Omega_x ; Omega_y ; Omega_z];  
        psi_z(i,j)= Rot(3,:)*[Omega_x ; Omega_y ; Omega_z];  
    end
end

end