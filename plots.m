%% Figures & Plots:
%% 2D plot of the distribution of Ct(k,t) (unaveraged)
% load("./results/19-Apr-2022_L_100_Jz_0.5_S_5_psi0_1_q0_0.31416.mat")
stag = 0;
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

for t = 1:RPT
%     ytotal = (L*norm(q)/psi_0)^2*(2*Ct{1}(1:L,t))/t^alpha;
%     ytotal =
%     (L*norm(q)/psi_0)^2*(2*C_t{1}(1:L,1,t)+C_t{3}(1:L,1,t))/T(t)^alpha;     %kx direction
    ytotal = (L*norm(q)/psi_0)^2*(2*Csum_t{1}(1,1:(L/2),t))/T(t)^alpha;
    xtotal = (2*pi*(0:(L/2-1))/L)/q(2)*T(t)^beta;
    plot(xtotal, ytotal,'Markersize',4,'Linewidth',1,'color',cc(t,:));
    hold on
end

set(gca,'XScale','log')
set(gca,'YScale','log')

if stag
    title(join(['$t^{',num2str(beta),'}k$-$\langle {S}_{-\bf k}^a(t){S}_{\bf k}^a(t)\rangle/t^{',num2str(alpha),'}\,q=',num2str(q(2)),'\,\mathrm{staggered\,spiral}$']))
else 
    title(join(['$t^{',num2str(beta),'}k$-$\langle {S}_{-\bf k}^a(t){S}_{\bf k}^a(t)\rangle/t^{',num2str(alpha),'}\,q=',num2str(q(2)),'\,\mathrm{spiral}$']))
end

%% 2D plot of the distribution of Ct(k,t): whole plane
load("concat/10-May-2022_L_400_Jz_0.5_S_10_psi0_0.70711_q0_0.31416_concat.mat")

surf(Csum_t{3}(:,:,3))

%%
%Plots the distribution of Ct(k,t): whole plane
axFtSz = 16 ; labFtSz = 16 ;
contourf(C_t{1}(:,:,90))
ax = gca;
ax.Box = 'on'
pbaspect([1,1,1])

%% Plots the distribution of Cav_t_xx(k,t): (radially averaged)
t_values = 1:RPT + 1;
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

T=[1:60];
s = size(T,2);
cc = jet(s);%cc: color legend 

for t = T
%     ytotal = (L*norm(q)/psi_0)^2*(2*Ct{1}(1:L,t))/t^alpha;
    ytotal = (L*norm(q)/psi_0)^2*(2*Cav_t{1}(1:L,t)+Cav_t{3}(1:L,t))/T(t)^alpha;
    xtotal = (2*pi*(0:L-1)/L)/q0*T(t)^beta;
    plot(xtotal([1:100]), ytotal([1:100]),'Markersize',4,'Linewidth',1,'color',cc(t,:));
    hold on
end

if stag
    title(join(['$t^{',num2str(beta),'}k$-$\langle {S}_{-\bf k}^a(t){S}_{\bf k}^a(t)\rangle/t^{',num2str(alpha),'}\,qN=',num2str(Nq),'\,\mathrm{staggered\,spiral}$']))
else 
    title(join(['$t^{',num2str(beta),'}k$-$\langle {S}_{-\bf k}^a(t){S}_{\bf k}^a(t)\rangle/t^{',num2str(alpha),'}\,qN=',num2str(Nq),'\,\mathrm{spiral}$']))
end


xlabel(['$t^\beta k/q_0$']);
ylabel(['$\langle {S}_{-\bf k}^a(t){S}_{\bf k}^a(t)\rangle/t^\alpha$']);
% set(gca,'XScale','log')
% set(gca,'YScale','log')
% xlim([0.01, 100]);
legend(int2str((1:RPT+1)'),'Location','southeast');
ax = gca;
ax.Box = 'on'
%% compute unequal time correlation

cc = jet(100);%cc: color legend 

t0 = 20
DT = 150

Ct_unequ = cell(1,3);
Ct_unequ(:) = {zeros(L,DT)};

C_t = cell(1,3); %individual components
C_t(:) = {zeros(L,L,DT)};%for plotting

Ct_k_test = cell(1,3); %individual components
Ct_k_test(:) = {zeros(L,DT)};%for plotting



axFtSz = 16 ; labFtSz = 16;
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultAxesFontName', 'Serif');
set(0,'defaultAxesFontSize',axFtSz);
set(0,'defaultTextFontSize',labFtSz);

for dt = 1 : DT
     for direc = 1 : 3 % for plotting
          C_t{direc}(:,:,dt) = abs(fft2(psi_t{direc}(:,:,t0))*conj(fft2(psi_t{direc}(:,:,dt+t0)))/L)^2; %2D cor relation function
          Ct_unequ{direc}(:,dt) = compute_C_unequ(psi_t{direc}(:,:,t0),psi_t{direc}(:,:,dt+t0),kr_range,nsigma,stag);
   
          Ct_k_test{direc}(:,dt) = compute_FT(psi_t{direc}(:,:,t0), psi_t{direc}(:,:,dt+t0));
     end
     disp('.')
%      ytotal = (L*norm(q)/psi_0)^2*(2*Ct_unequ{1}(1:L,dt) + Ct_unequ{3}(1:L,dt));
%      ytotal = abs((L*norm(q)/psi_0)^2*(2*Ct_k_test{1}(1:L,dt) + Ct_k_test{3}(1:L,dt)));
%      xtotal = (2*pi*(0:L-1)/L);
%      plot(xtotal([1:100]), ytotal([1:100]),'Markersize',4,'Linewidth',1,'color',cc(dt,:));
%      hold on
end

%% plot fourier transformed unequal time correlation

% wrange = (0:(2*size(Ct_k_test{direc}(:,1:DT),2)-2))*2*pi/(size(Ct_k_test{direc}(:,1:DT),1)*tau);
% [C,h] =contourf(1:L,wrange,transpose(log(abs(compute_fft(Ct_k_test{direc}(:,1:DT),tau)))),40);
[C,h] = contourf(transpose(log(abs(compute_fft(Ct_k_test{direc}(:,1:DT),tau)))),40);
set(h,'LineColor','none');
xlabel(['$k$']);
ylabel(['$\omega\tau_\ast$']);
pbaspect([1,1,1]); 
c = colorbar('location','northoutside')
c.Label.String = '$\langle {S}_{-\bf k}^a(t){S}_{\bf k}^a(t)\rangle$';
c.Label.Interpreter = 'latex'
% set(gca,'XScale','log')
%%
contourf(transpose(log(abs(Ct_k_test{direc}))))
%% plot unequal time correlation


cc = jet(RPT - t0);%cc: color legend 

axFtSz = 16 ; labFtSz = 16 ;
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultAxesFontName', 'Serif');
set(0,'defaultAxesFontSize',axFtSz);
set(0,'defaultTextFontSize',labFtSz);

for dt =  1 : 90
%     ytotal = (L*norm(q)/psi_0)^2*(2*Ct{1}(1:L,t))/t^alpha;
    ytotal = (L*norm(q)/psi_0)^2*(2*Ct_unequ{1}(1:L,dt) + Ct_unequ{3}(1:L,dt));
    xtotal = (2*pi*(0:L-1)/L);
    plot(xtotal([1:100]), ytotal([1:100]),'Markersize',4,'Linewidth',1,'color', cc(dt,:));
    hold on
end

if stag
    title(join(['$t^{',num2str(beta),'}k$-$\langle {S}_{-\bf k}^a(t){S}_{\bf k}^a(t)\rangle/t^{',num2str(alpha),'}\,qN=',num2str(Nq),'\,\mathrm{staggered\,spiral}$']))
else 
    title(join(['$t^{',num2str(beta),'}k$-$\langle {S}_{-\bf k}^a(t){S}_{\bf k}^a(t)\rangle/t^{',num2str(alpha),'}\,qN=',num2str(Nq),'\,\mathrm{spiral}$']))
end

xlabel(['$t^\beta k/q_0$']);
ylabel(['$\langle {S}_{-\bf k}^a(t_0){S}_{\bf k}^a(t)\rangle/t^\alpha$']);
% set(gca,'XScale','log')
% set(gca,'YScale','log')
legend(int2str((1:RPT + 1 - t0)'),'Location','southeast');
ax = gca;
ax.Box = 'on'
xlim([0.1, 10]);
ylim([10^2, 8E5]);

%%
% contourf(log(transpose(abs(Ct_unequ{direc}))),30)
contourf(log(transpose(abs(fft(Ct_unequ{direc})))),30)
% set(gca,'XScale','log')
% values = hist3([data1(:) data2(:)],[51 51]);
% imagesc(values.')

%%
%Plots spin vectors on a two-dimensional plane for different times
h = figure;

axFtSz = 16 ; labFtSz = 16 ;
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultAxesFontName', 'Serif');
set(0,'defaultAxesFontSize',axFtSz);
set(0,'defaultTextFontSize',labFtSz) ;
ax = gca;
ax.Box = 'on';
set(gca,'linewidth',3,'XColor',[0 0 0],'YColor',[0 0 0])

axis tight manual % this ensures that getframe() returns a consistent size
filename = 'quiverAnimated.gif';


for i = 1
    quiver(1:L,1:L,psi_t{1}(:,:,i),psi_t{2}(:,:,i),'Linewidth',1,'ShowArrowHead','On')
    pbaspect([1,1,1])
    axis([-6 L + 6 -6 L+6])
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

%%
%Plots spin density on a two-dimensional plane for different times
axFtSz = 16 ; labFtSz = 16 ;
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultAxesFontName', 'Serif');
set(0,'defaultAxesFontSize',axFtSz);
set(0,'defaultTextFontSize',labFtSz) ;

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'ContourAnimated.gif';
for i = 1:30
    [C,y]= contourf(psi_t{3}(:,:,i));
     caxis([-1 1])
     colorbar
    set(y,'LineColor','none');
    pbaspect([1,1,1])
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

%%
%Plots Ct_xx density on a two-dimensional plane:
axFtSz = 16 ; labFtSz = 16 ;
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultAxesFontName', 'Serif');
set(0,'defaultAxesFontSize',axFtSz);
set(0,'defaultTextFontSize',labFtSz) ;

h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = 'Ct_ii_Animated.gif';

for i = 1:RPT+1
    % Draw plot for y = x.^n
    surf((2*pi*(1:L)/L),(2*pi*(1:L)/L),C_t{1}(:,:,i))
    zlim([0 200])
%     axis([-1 2*pi+1 -1 2*pi+1 ])
    pbaspect([1,1,1])
    caxis([0, 250]);
    colorbar;
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if i == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
end

%%
% Fitting scaling function
%Define scaling parameters 
R = 200;
Ra = 30;
Rb = 30;
e = zeros(Ra,Rb);
amin=0;
amax=2;
bmin=0;
bmax=1;
valpha=[amin:(amax-amin)/(Ra):amax];
vbeta=[bmin:(bmax-bmin)/(Rb):bmax];
%sample range

for i = 1:8
    t0 = 10 + i*5
    t1 = t0 - 10;
    t2 = t0 + 10;
    for a = 1:Ra
        for b = 1:Rb
            alpha = valpha(a);
            beta = vbeta(b);
            xtotal = ((2*pi*(0:L-1)/L)/q0)';
            x = xtotal([1:R/2]);


            for t =  t1 : t2 - 1
                for s =  t1 : t2 - 1     
%                    ytotal1 = (L*norm(q)/psi_0)^2*(2*Ct{1}(1:R,t)+Ct{3}(1:R,t));
%                    ytotal2 = (L*norm(q)/psi_0)^2*(2*Ct{1}(1:R,s)+Ct{3}(1:R,s));

                    %For Joaquin's Data
                   ytotal1 = 2*Cav_t{1}(1:R,t)+Cav_t{3}(1:R,t);
                   ytotal2 = 2*Cav_t{1}(1:R,s)+Cav_t{3}(1:R,s);

                   y1 = ytotal1([1:R/2]);
                   y2 = ytotal2([1:R/2]);

                   xq1 = xtotal([1:R/2])/(t/t0)^beta;
                   xq2 = xtotal([1:R/2])/(s/t0)^beta;

                   yq1 = interp1(x,y1,xq1,'spline')/(t/t0)^alpha; %standard
                   yq2 = interp1(x,y2,xq2,'spline')/(s/t0)^alpha;

                     dy = abs(yq1-yq2);
                     Sq = sum(dy); 
                     e(a,b) = e(a,b) + Sq;

                end
            end  
        end
    end
    [maxValue, linearIndexesOfMaxes] = min(e(:))
    mod(linearIndexesOfMaxes,Ra)/Ra*(amax-amin)+amin
    floor(linearIndexesOfMaxes/Ra)/Rb*(bmax-bmin)+bmin
    
    subplot(1, 8 ,i)
    contourf(vbeta([1:Rb]),valpha([1:Ra]),log(e));
    pbaspect([1 1 1])
    xlabel(['$\beta$']);
    if i == 1
         ylabel(['$\alpha$']);
    end
    if i > 1     
        set(gca,'ytick',[])
    end
    hold on
end
% title(join('Error function $t_0=$',num2str(t0)));


%% Exporting Data
to_save = (L*norm(q)/psi_0)^2*(2*Cav_t{1}(:,:,:)+Cav_t{3}(:,:,:))
save('M.mat','to_save')

%%