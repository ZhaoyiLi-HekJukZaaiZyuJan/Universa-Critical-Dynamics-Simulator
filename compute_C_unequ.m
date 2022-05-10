function C = compute_C_unequ(psi_i_t0, psi_i,kr_range,nsigma,stag)
L = size(psi_i,1);

%Compute FT:
FT = fft2(psi_i)/L;
FT_t0 = fft2(psi_i_t0)/L;


%Average over wavevectors:
k_range = -pi+2*pi/L:2*pi/L:pi;
sigma = nsigma*2*pi/L;
C = zeros(L,1);

%1:


for r = 1 : L
    if stag == 0
        kr = k_range(mod2(r+ceil(L/2)-1,L)); %r=1 k=0/r=51, k=pi/r=L k=-0.0628
    end
    if stag == 1
        kr = k_range(mod2(r-1,L));
    end
    
    
    norm = 0;
    for i = 1 : L
        ki = k_range(mod2(i+ceil(L/2)-1,L));
        for j = 1 : L
            kj = k_range(mod2(j+ceil(L/2)-1,L));
            if stag == 0
                weight = exp(-(sqrt(ki^2+kj^2)-abs(kr))^2/sigma^2);
            end
            if stag == 1
                weight = exp(-(sqrt((abs(ki)-pi)^2+(abs(kj)-pi)^2)-abs(abs(kr)-pi))^2/sigma^2);
            end
            C(r) = C(r) + real(FT(i,j)*conj(FT_t0(i,j)))*weight;
            norm = norm + weight;
        end
    end
    C(r) = C(r)/norm;
    
end

end
