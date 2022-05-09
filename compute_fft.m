function Fw_xx = compute_fft(Ft_xx,tau)

%Compute FFT with frequencies on the scale of J:

k_values = size(Ft_xx,1);
T = size(Ft_xx,2);
Ts = 2*T - 1;
y = zeros(1,Ts+1);
Fw_xx = zeros(k_values,Ts+2);

for k = 1 : k_values
    y(T+1:Ts+1) = (Ft_xx(k,:));
    y(1:T) = 0*(flip(Ft_xx(k,:)));
    Y_ft0 = fft([0,y]);
    Fw_xx(k,:) = Y_ft0;
end
omega_range = 0:1:(Ts-1);
omega_range = omega_range*2*pi()/(Ts*tau);

end

