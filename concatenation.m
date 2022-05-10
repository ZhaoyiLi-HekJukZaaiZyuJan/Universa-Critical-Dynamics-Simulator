%%concatenate all simulated data
name = "random_24-Apr-2022_L=100_Jz=0.5_S=10_psi0=0.70711_stag=1_q0=0.31416";

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

filename = strcat('/scratch/users/ladmon/AFM/results/',name, '_concat.mat');
save(filename,'Csum_t','L','J_z','S','psi_0','q','q0','Nq','stag','N_samples','kr_range','T','tau');
%save(filename,'Csum_t','L','J_z','S','psi_0','q','N_samples','kr_range','T','tau');
