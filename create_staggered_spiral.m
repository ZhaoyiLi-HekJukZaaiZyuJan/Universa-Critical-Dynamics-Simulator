<<<<<<< HEAD
function [psi_x,psi_y,psi_z] = create_staggered_spiral(L,qx,qy,t,S)
=======
function [psi_x,psi_y,psi_z] = create_staggered_spiral(L,qx,qy,t,S,stag)
>>>>>>> origin/master

psi_x = zeros(L,L);
psi_y = zeros(L,L);
psi_z = zeros(L,L);

for j = 1 : L
    for k = 1 : L
            phi = qx*j+ qy*k;
<<<<<<< HEAD
            theta = t * (-1)^(j+k);
            if mod(j+k,2)
=======

            if mod(j+k,2) | stag == 0
>>>>>>> origin/master
               theta = t;
            else
               theta = t+pi;
            end

        
            Rot = [cos(theta)*cos(phi) -sin(phi) sin(theta)*cos(phi); cos(theta)*sin(phi) cos(phi) sin(theta)*sin(phi);-sin(phi) 0 cos(theta)];
            S_x = normrnd(0, 1/sqrt(2*S));
            S_y = normrnd(0, 1/sqrt(2*S));
%             S_z = normrnd(S, 1/(2*sqrt(2)));
            S_z = 1;
%             
%             psi_x(j,k) = Rot(1,:)*[S_x;S_y;S_z];  
%             psi_y(j,k) = Rot(2,:)*[S_x;S_y;S_z];  
%             psi_z(j,k) = Rot(3,:)*[S_x;S_y;S_z];  

            
            psi_x(j,k) = Rot(1,:)*[S_x;S_y;S_z];  
            psi_y(j,k) = Rot(2,:)*[S_x;S_y;S_z];  
            psi_z(j,k) = Rot(3,:)*[S_x;S_y;S_z];  

    end
end

end