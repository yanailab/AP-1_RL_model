% Loop over realizations and parameters
%close all; clear all;

Ms   = [5 10 15 20]
Ns   = [12 14 16 18 20]
% 
% Ms   = [5 10]
% Ns   = [12 16 20]

%stats = zeros(length(Ms),10);
S_stats = zeros(length(Ms),length(Ns));
T_stats = zeros(length(Ms),length(Ns));
M_stats = zeros(length(Ms),length(Ns));

    Nloops = 50;   % Number of stochastic realizations to average over


for kM = 1:length(Ms)

    for kN = 1:length(Ns)
  
    M = Ms(kM);
    N = Ns(kN);


Model3_3state_DiscreteM_ForLoops;  %Run once to set sizes

Stress_loop = zeros(length(time),Nloops);
B_loop_plus = zeros(length(time),Nloops);
B_loop_mnus = zeros(length(time),Nloops);
Pres_loop   = zeros(length(time),Nloops);
CondAc_loop = zeros(length(time),Nloops);

Stress_loop(:,1)      = res.stress;
B_loop_plus(:,1)      = mBp;
B_loop_mnus(:,1)      = mBm;
Pres_loop(:,1)        = mP;
CondAc_loop(:,1)      = frac_active_given_present;

for jj=2:Nloops
    rng('shuffle');
    Model3_3state_DiscreteM_ForLoops;
    Stress_loop(:,jj) = res.stress;
    B_loop_plus(:,jj)      = mBp;
    B_loop_mnus(:,jj)      = mBm;
    Pres_loop(:,jj)        = mP;
    CondAc_loop(:,jj)      = frac_active_given_present;
    %plot(B_keep(:,1:N/2), 'LineWidth',1,'Color',[0 0.4470 0.7410]); hold on
    %plot(B_keep(:,N/2+1:N), 'LineWidth',1,'Color',[0.85 0.325 0.098]);

end

figure;   % STRESS NORMALIZED BY f
sig = std(Stress_loop');
errorbar(time',mean(Stress_loop,2)./f,sig./f,'Color',[0 0.4470 0.7410]*1.2);
hold on;
plot(time,mean(Stress_loop,2)./f,'-k','LineWidth',2);
xlabel('Time Steps'); ylabel('Stress');
title(num2str(M));

ms = mean(Stress_loop,2)./f;
S_stats(kM,kN) = mean(ms(floor(end*0.75):end));
idx = find(smooth(ms)-0.5>0);
half = idx(end);
t_half = time(half);
%PP=polyfit(time(1:t_fit),log(ms(1:t_fit)),1);
T_stats(kM,kN) = t_half;
% figure;plot(time,mean(Stress_loop,2)./f-S_stats(kM,kN));
% hold on;
% plot(time,PP(1)*time+PP(2));


% figure
% sigp = std(B_loop_plus');
% errorbar(time',mean(B_loop_plus,2),sigp,'Color',[0 0.4470 0.7410]*1.2);
% hold on;
% sigm = std(B_loop_mnus');
% errorbar(time',mean(B_loop_mnus,2),sigp,'Color',[0.85 0.325 0.098]);
 mb1 = mean(B_loop_plus,2);
 M_stats(kM,kN) = mean(mb1(floor(end/2):end))/M;
%stats(kM,3) = mean(mb1(floor(end/2):end))/M;
%stats(kM,4) = std(mb1(floor(end/2):end))/M;
%[mean(mb1(floor(end/2):end)) std(mb1(floor(end/2):end))]
%mb2 = mean(B_loop_mnus,2);
%stats(kM,5)=mean(mb2(floor(end/2):end))/M;
%stats(kM,6)=std(mb2(floor(end/2):end))/M;
%[mean(mb2(floor(end/2):end)) std(mb2(floor(end/2):end))]
%title(num2str(M));

% figure;
% sig = std(Pres_loop');
% errorbar(time',mean(Pres_loop,2),sig,'Color',[0 0.4470 0.7410]*1.2);
% hold on;
% sig = std(CondAc_loop');
% plot(time,mean(CondAc_loop,2),'-k','LineWidth',2);
% xlabel('Time Steps'); legend ('Fraction Present','Active | Present')
% title(num2str(M));

%p1 = mean(Pres_loop,2);
%stats(kM,7)=mean(p1(floor(end/2):end)) ;
%stats(kM,8)=std(p1(floor(end/2):end));
%[mean(p1(floor(end/2):end)) std(p1(floor(end/2):end))]
%p2 = mean(CondAc_loop,2);
%stats(kM,9) = mean(p2(floor(end/2):end));
%stats(kM,10) = std(p2(floor(end/2):end));
%[mean(p2(floor(end/2):end)) std(p2(floor(end/2):end))]

    end
end

%%%%%%%%
% PLOTS
%%%%%%%%

figure;imagesc(Ns,Ms,S_stats);axis xy; axis square; colormap(bone);
shading flat;
%clim([0 0.12]); 
colorbar
xlabel('$N$','Interpreter', 'latex');ylabel('$M$','Interpreter', 'latex');
ax = gca;                     
ax.FontSize = 14; 

figure;imagesc(Ns,Ms,T_stats);axis xy; axis square; colormap(bone);
shading flat;
%clim([0 0.12]); 
colorbar
xlabel('$N$','Interpreter', 'latex');ylabel('$M$','Interpreter', 'latex');
ax = gca;                     
ax.FontSize = 14; 
%exportgraphics(gca,'Stress_Loop10_Time40K.png','Resolution',600)


figure;imagesc(Ns,Ms,M_stats);axis xy; axis square; 
%colormap(flipud(pink));
colormap(bone);
shading flat;
clim([0.2 1]); 
colorbar
xlabel('$N$','Interpreter', 'latex');ylabel('$M$','Interpreter', 'latex');
ax = gca;                     
ax.FontSize = 14; 
%exportgraphics(gca,'Specificity_Loop10_Time40K.png','Resolution',600)


% figure;plot(Ns,S_stats,'-*','linewidth',2);xlabel('N');ylabel('Reduced Stress');
% figure;plot(Ms,S_stats','-*','linewidth',2);xlabel('M');ylabel('Reduced Stress');
% 
% figure;plot(Ns,M_stats,'-*','linewidth',2);xlabel('N');ylabel('Stability');
% figure;plot(Ms,M_stats','-*','linewidth',2);xlabel('M');ylabel('Stability');
 
M_all = M_stats(:);
S_all = S_stats(:);
figure;
plot(M_all,log((1-S_all)./S_all),'ok', 'linewidth',2,'MarkerFaceColor',[0.5,0.5,0.5]);
xlabel('$m^*/M$','Interpreter', 'latex');ylabel('$\ln (1-S^*)/S^*$','Interpreter', 'latex');
ax = gca;                     
ax.FontSize = 14; 
%exportgraphics(gca,'M_S_Loop20_Time40k.png','Resolution',600)




% figure;
% for ii=1:length(Ns)
%     plot(Ns(ii),T_stats(2,ii),'v','LineWidth',3);hold on
% end;
% 
% figure;
% for jj=1:length(Ms)
%     plot(Ms(jj),T_stats(jj,1),'v','LineWidth',3);hold on
%     plot(Ms(jj),T_stats(jj,2),'o','LineWidth',3)
