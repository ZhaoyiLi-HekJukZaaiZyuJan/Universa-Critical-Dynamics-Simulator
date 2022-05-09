%%concatenate all simulated data
name = "21-Apr-2022_L=100_Jz=0.5_S=10_psi0=0.70711_q0=0.62832";

load(strcat("./results/",name,"_i=1.mat"));
RPT = size(T,2);

Csum_t = cell(1,3);
Csum_t(:) = {zeros(L,L,RPT+1)};

a=0;
for i = [1:100]
    disp('.')
    try
        load(strcat("./results/", name, "_i=",num2str(i),".mat"))
        for direc = [1 2 3]
            Csum_t{direc}(:,:,:) = Csum_t{direc}(:,:,:) + C_t{direc}(:,:,:);%2D correlation function
        end
        a = a+1;
    end

end

filename = strcat('./results/',name, '_concat.mat');
save(filename,'Csum_t','L','J_z','S','psi_0','q','N_samples','kr_range','T','tau');

%% 2D plot of the distribution of Ct(k,t) (unaveraged)

%load concatenated%
%load(strcat("./results/21-Apr-2022_L=100_Jz=0.5_S=10_psi0=0.70711_stag=1_q0=0.31416_concat.mat"));
RPT = size(T,2);
t_values = 1:RPT;
% k_inertial = kr_range(kmin:kmax)

%Define time values

%Define scaling parameters 
% 
alpha = 0;
beta = 0;
%theta = asin(psi_0)

axFtSz = 16 ; labFtSz = 16 ;
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultAxesFontName', 'Serif');
set(0,'defaultAxesFontSize',axFtSz);
set(0,'defaultTextFontSize',labFtSz);

cc = jet(RPT);%cc: color legend 

for t = 2:100
    ytotal = (L*norm(q)/psi_0)^2*(2*Csum_t{2}(1+stag*(L/2),(1:L/2)+stag*(L/2),t))/T(t)^alpha;
    xtotal = (2*pi/L*(0:(L/2-1)))/q(2)*T(t)^beta;
    plot(xtotal, ytotal,'Markersize',4,'Linewidth',1,'color',cc(t,:));
    hold on
end
% 
% set(gca,'XScale','log')
% set(gca,'YScale','log')

if stag
    title(join(['$t^{',num2str(beta),'}k$-$\langle {S}_{-\bf k}^a(t){S}_{\bf k}^a(t)\rangle/t^{',num2str(alpha),'}\,qN=',num2str(Nq),'\,\mathrm{staggered\,spiral}$']))
else 
    title(join(['$t^{',num2str(beta),'}k$-$\langle {S}_{-\bf k}^a(t){S}_{\bf k}^a(t)\rangle/t^{',num2str(alpha),'}\,qN=',num2str(Nq),'\,\mathrm{spiral}$']))
end
%%
% 2D plot of the distribution of Ct(k,t): whole plane
% load("results/19-Apr-2022_L_100_Jz_0.5_S_10_psi0_0.70711_q0_0.62832.mat")
fname = "random_24-Apr-2022_L=100_Jz=0.5_S=10_psi0=0.70711_stag=1_q0=0.31416";
%fname = "random_24-Apr-2022_L=100_Jz=0.5_S=10_psi0=0.70711_stag=1_q0=0.31416",
load(strcat("./results/", fname,"_concat.mat"))
% surf(Csum_t{1}(:,:,1))
