function [psi_xout,psi_yout,psi_zout] = run_twark4(psi_x,psi_y,psi_z,J_z,T,N_steps)

t_range = 0 : T/N_steps : T;
dt = T/N_steps;
L = size(psi_x,1);
psi_x2 = zeros(L,L);
psi_y2 = zeros(L,L);
psi_z2 = zeros(L,L);
Hx = zeros(L,L);
Hy = zeros(L,L);
Hz = zeros(L,L);


psi_xout = psi_x;
psi_yout = psi_y;
psi_zout = psi_z;
psi_x0 = psi_x;
psi_y0 = psi_y;
psi_z0 = psi_z;
for k = 1 : N_steps
    psi_x2 = psi_x0;
    psi_y2 = psi_y0;
    psi_z2 = psi_z0;
    [Hx,Hy,Hz] = Derivative(psi_x2,psi_y2,psi_z2,J_z);
    psi_xout = psi_xout + dt*Hx/6;
    psi_yout = psi_yout + dt*Hy/6;
    psi_zout = psi_zout + dt*Hz/6;

    psi_x2 = psi_x0 + dt*Hx/2;
    psi_y2 = psi_y0 + dt*Hy/2;
    psi_z2 = psi_z0 + dt*Hz/2;
    [Hx,Hy,Hz] = Derivative(psi_x2,psi_y2,psi_z2,J_z);
    psi_xout = psi_xout + dt*Hx/3;
    psi_yout = psi_yout + dt*Hy/3;
    psi_zout = psi_zout + dt*Hz/3;

    psi_x2 = psi_x0 + dt*Hx/2;
    psi_y2 = psi_y0 + dt*Hy/2;
    psi_z2 = psi_z0 + dt*Hz/2;
    [Hx,Hy,Hz] = Derivative(psi_x2,psi_y2,psi_z2,J_z);
    psi_xout = psi_xout + dt*Hx/3;
    psi_yout = psi_yout + dt*Hy/3;
    psi_zout = psi_zout + dt*Hz/3;

    psi_x2 = psi_x0 + dt*Hx;
    psi_y2 = psi_y0 + dt*Hy;
    psi_z2 = psi_z0 + dt*Hz;
    [Hx,Hy,Hz] = Derivative(psi_x2,psi_y2,psi_z2,J_z);
    psi_xout = psi_xout + dt*Hx/6;
    psi_yout = psi_yout + dt*Hy/6;
    psi_zout = psi_zout + dt*Hz/6;
    
    psi_x0 = psi_xout;
    psi_y0 = psi_yout;
    psi_z0 = psi_zout;    
end

end