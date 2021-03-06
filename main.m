%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main.m %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Created by Joaquin, Modified by Zhaoyi
%% Part 1. Simulation:
function main(seed)
%%
seed = 1 %for local simulation
%input parameters:
L = 100;
kr_range = 2*pi/L:2*pi/L:pi;
S = 1;
theta = pi/3;
nsigma = 2; % gaussian initial condition

%Derived parameters:

N_samples = 1; 
J_z = 0.5;
Nq = 20;
q0 = Nq*2*pi/L;
q = [0, q0];
psi_0 = sin(theta);
stag = 0; %1 (anti) or 0 (ferro)
random = 1;
% t0 = 20; %unequ time



% Time range:
tau = 1/(psi_0^2*(2-cos(q(1))-cos(q(2))));%defines timescale of the spiral state
% Stepsize：
% T = [tau*ones(1,20),tau/15*ones(1,150)];
T = [tau*ones(1,2)];%Stepsize
% T = [tau*ones(1,60)];%Stepsize
RPT = size(T,2);
N_steps = 100;%Roughly speaking, you want to do ~1000-2000steps/tau
%initialize:


Magnetization = zeros(RPT+1,1);
Energy = zeros(RPT+1,1);
% Cav_t = cell(1,3);
% Cav_t(:) = {zeros(L,RPT+1)};

%non-averaging
C_t = cell(1,3);
C_t(:) = {zeros(L,L,RPT+1)};


fprintf('---Start Sampling---\n\n')
for sample = 1 : N_samples
    %Create sample:
    psi_t = cell(1,3);
    psi_t(:) = {zeros(L,L,RPT+1)}; %for plotting, comment out to save memory
    psi = cell(1,3); %at particular times
    psi(:) = {zeros(L,L)};
    
    if random == 1
        [psi{1},psi{2},psi{3}] = create_random(L,theta,S,seed,stag);
    else
          [psi{1},psi{2},psi{3}] = create_staggered_spiral(L,q(1),q(2),theta,S,stag);
    end
    
    
    %run
    for i = 1 : RPT
        if i ~= 1
            tic
            [psi{1},psi{2},psi{3}] = run_twark4(psi{1},psi{2},psi{3},J_z,T(i), N_steps);%step foward
        end
        
        for direc = 1 : 3 % for plotting
            psi_t{direc}(:,:,i) = psi{direc};
        end
        
        % Compute evolved quantities:
        Magnetization(i) = Magnetization(i) + compute_magnetization(psi{3});
        Energy(i) = Energy(i) + compute_energy(psi{1},psi{2},psi{3},J_z);
        
        % Compute evolved correlations:
        for direc = [1 2 3]
%           Cav_t{direc}(:,i) = Cav_t{direc}(:,i) + compute_C_direct(psi{direc},kr_range,nsigma,stag);
            %Ct_unequ{direc}(:,i) = Ct_unequ{direc}(:,i) + C_unequ; 
            C_t{direc}(:,:,i) = C_t{direc}(:,:,i) + abs(fft2(psi{direc})/L).^2;%2D correlation function
        end
        if i ~= 1
            toc
        end
        fprintf('\nSampling percentage: %6.2f.\n',sample/N_samples*100)
        fprintf('Step %6.2f.\n',i/(RPT+1)*100)
    end
end

%Normalization:
Energy = Energy/N_samples;
Magnetization = Magnetization/N_samples;
for direc = [1 2 3]
%     Cav_t{direc} =  Cav_t{direc}/N_samples;
    C_t{direc} = C_t{direc}/N_samples; %change this
end
%get date for filename
c=date;

filename = strcat('/scratch/users/ladmon/AFM/results/',c, '_L_',num2str(L), '_Jz_', num2str(J_z), ...
         '_S_',num2str(S),'_psi0_', num2str(psi_0), '_q0_', num2str(q0),'_i=',num2str(seed),'.mat');

% save(filename,'psi_t','C_t','Energy','Magnetization','L','J_z','S','psi_0','q','q0','Nq','N_samples','kr_range','stag','T','tau');
disp("saved")
%% 2Splot C_t, check

surf(C_t{3}(:,:,1))
