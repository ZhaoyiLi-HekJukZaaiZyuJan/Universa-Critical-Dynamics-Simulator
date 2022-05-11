%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Concatenation.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%,

%%concatenate all simulated data
name = "10-May-2022_L_400_Jz_0.5_S_10_psi0_0.70711_q0_0.31416";

load(strcat("/scratch/users/ladmon/AFM//results/",name,"_i=1.mat"));

RPT = size(T,2);

Csum_t = cell(1,3);
Csum_t(:) = {zeros(L,L,RPT+1)};

a=0;
for i = [1:100]
    disp('.')
    try
        load(strcat("/scratch/users/ladmon/AFM/results/", name, "_i=",num2str(i),".mat"))
        for direc = [1 2 3]
            Csum_t{direc}(:,:,:) = Csum_t{direc}(:,:,:) + C_t{direc}(:,:,:);%2D correlation function
        end
        a = a+1;
    end

end

filename = strcat('/scratch/users/ladmon/AFM/results/concat/',name, '_concat.mat');
filename = strcat('~/AFM/concat/',name, '_concat.mat');
save(filename,'Csum_t','L','J_z','S','psi_0','N_samples','kr_range','T','tau');
%save(filename,'Csum_t','L','J_z','S','psi_0','q','N_samples','kr_range','T','tau');
