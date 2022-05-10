%Function that creates a spiral state: the ansats I use is that, if
%S_x^2+S_y^2>1, then I allow the S_z component to be zeros using the function
%f(x)=x if x<1, and f(x)=2-x if 1<x<2, and f(x) = 0 if x>2.

<<<<<<< HEAD
function [psi_x,psi_y,psi_z] = create_spiralgaussian(L,qx,qy,theta,S,seed)
=======
function [psi_x,psi_y,psi_z] = create_spiralgaussian(L,qx,qy,t,S,seed,stag)
>>>>>>> origin/master

	psi_x = zeros(L,L);
	psi_y = zeros(L,L);
	psi_z = zeros(L,L);
	
<<<<<<< HEAD
=======
            
>>>>>>> origin/master
	for i = 1 : L
		for j = 1 : L
			Omega_x = normrnd(0, 1/sqrt(2*S));
			Omega_y = normrnd(0, 1/sqrt(2*S));
			Omega_z = 1;
<<<<<<< HEAD
			disp(Omega_x);
=======
>>>>>>> origin/master
			psi_x(i,j)=Omega_x;
			psi_y(i,j)=Omega_y;
			psi_z(i,j)=Omega_z;
		end
<<<<<<< HEAD
	end
    rng(1)
	for i = 1 : L
		for j = 1 : L
			
			phi =  2*pi*rand();
			Rot = [cos(theta)*cos(phi) -sin(phi) sin(theta)*cos(phi);cos(theta)*sin(phi) cos(phi) sin(theta)*sin(phi);-sin(phi) 0 cos(theta)];
            Omega_x = psi_x(i,j);
            Omega_y = psi_y(i,j);
            Omega_z = psi_z(i,j);
=======
     end
     
	for i = 1 : L
		for j = 1 : L
			if mod(i+j,2) | stag == 0
                    theta = t;
               else
                    theta = t+pi;
               end
      
			phi =  2*pi*rand();
			Rot = [cos(theta)*cos(phi) -sin(phi) sin(theta)*cos(phi);cos(theta)*sin(phi) cos(phi) sin(theta)*sin(phi);-sin(phi) 0 cos(theta)];
               Omega_x = psi_x(i,j);
               Omega_y = psi_y(i,j);
               Omega_z = psi_z(i,j);
>>>>>>> origin/master
			psi_x(i,j)= Rot(1,:)*[Omega_x ; Omega_y ; Omega_z];
			psi_y(i,j)= Rot(2,:)*[Omega_x ; Omega_y ; Omega_z];
			psi_z(i,j)= Rot(3,:)*[Omega_x ; Omega_y ; Omega_z];

		end
	end

			
	
	end