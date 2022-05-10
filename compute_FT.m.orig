function Sxx = compute_FT(psi_x0,psi_x)

L = size(psi_x,1);

A = fft2(psi_x0);
B = fft2(psi_x);

shift = round(L/2)-1;
FT = zeros(L,L);
Sxx = zeros(L,1);
Mxx = zeros(L,1);

for i = 1 : L
    i_new = mod2(i+shift,L);
    for j = 1 : L
        j_new = mod2(j+shift,L);
        FT(i_new,j_new) = conj(A(i,j))*B(i,j);
    end
end

kx_range = (-pi+2*pi()/L):(2*pi()/L):(pi());

for i = 1 : L
    n_y = round(L/2);
    Sxx(i) = FT(i,n_y);
    n_y = L;
    Mxx(i) = FT(i,n_y);
end


end