%Function that creates a spiral state: the ansats I use is that, if
%S_x^2+S_y>1, then I allow the S_z component to be zeros using the function
%f(x)=x if x<1, and f(x)=2-x if 1<x<2, and f(x) = 0 if x>2.

function [psi_x,psi_y,psi_z] = create_spiralgaussian(L,qx,qy,theta,S)

psi_x = zeros(L,L);
psi_y = zeros(L,L);
psi_z = zeros(L,L);

for i = 1 : L
    for j = 1 : L
        Omega_x = normrnd(0,1/sqrt(2*S));
        Omega_y = normrnd(0,1/sqrt(2*S));
        disp(Omega_x);
        Omega_z = 1;

        phi = qx*i+qy*j;
        Rot = [cos(theta)*cos(phi) -sin(phi) sin(theta)*cos(phi);cos(theta)*sin(phi) cos(phi) sin(theta)*sin(phi);-sin(phi) 0 cos(theta)];
        psi_x(i,j)= Rot(1,:)*[Omega_x ; Omega_y ; Omega_z];  
        psi_y(i,j)= Rot(2,:)*[Omega_x ; Omega_y ; Omega_z];  
        psi_z(i,j)= Rot(3,:)*[Omega_x ; Omega_y ; Omega_z];  
    end
end

end