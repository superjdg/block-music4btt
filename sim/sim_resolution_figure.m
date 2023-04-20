%20230420 plot resolution figure
%Fig 9 in this article. 
%@article{wang2023block, 
% title={Block-MUSIC in Blade Tip Timing: Performance Study of Block Snapshot Matrix}, 
% author={Wang, Zengkun and Yang, Zhibo and Wu, Shuming and Tian, Shaohua and Chen, Xuefeng}, 
% journal={Mechanical Systems and Signal Processing}, year={2023}}
clear all
close all
load Prob_G; 
%% plot
colors = ["#1A2B3E","#185F8C","#97BEA1","#A4A750","#8B5C0C"];
linestyles = ["-","--","-.",":"];
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
nexttile
hold on
for i = 1:4
    plot(normalized_f,Prob_G(i,:),linestyles(i),'Color',colors(i));
end
xlim([min(normalized_f),max(normalized_f)]);
ylim([0, 1]);
box on
hl = legend("$G=1$","$G=2$","$G=4$",...
    "$G=8$","$G=16$",'Interpreter','latex','location','north','numcolumns',1);
set(hl,'FontName','Times New Roman','FontSize',6,'FontWeight','normal')
set(hl,'Box','off');
xlabel("$\omega_r-\omega_l$",'Interpreter','latex')
ylabel("Probability");
set(gca,'FontName','Times New Roman')
title('(a)');

load P_G_Rayleighs.mat;
nexttile
hold on
semilogy(n_ps,P_rec_G,...
    linestyles(1),'Color',colors(1));
hold on
semilogy(n_ps,Rayleighs, ...
    linestyles(2),'Color',colors(2));
xlim([min(n_ps),max(n_ps)]);
% ylim([0,max(Rayleighs)]);
xlabel("$G$",'Interpreter','latex')
ylabel("$\omega_r-\omega_l$",'Interpreter','latex');
set(gca,'FontName','Times New Roman')
ax = gca;
ax.YAxis.Exponent = -3;
hl = legend(['$M=80$',newline,'$Q=34$'],['Rayleigh',newline,'length'],...
    'Interpreter','latex','location','south','numcolumns',1);
set(hl,'Box','off');
box on

set(gcf,'unit','centimeters','position',[3 5 8 6])
title('(b)');



    
    