%20230420 
%plot Fig 7 in this article. 
%@article{wang2023block, 
% title={Block-MUSIC in Blade Tip Timing: Performance Study of Block Snapshot Matrix}, 
% author={Wang, Zengkun and Yang, Zhibo and Wu, Shuming and Tian, Shaohua and Chen, Xuefeng}, 
% journal={Mechanical Systems and Signal Processing}, year={2023}}
clear all
close all

load E_amax_amin.mat; load E_G; 
load E_M ;load E_sigma_amin;


%% »­Í¼
colors = ["#1A2B3E","#185F8C","#97BEA1","#A4A750","#8B5C0C"];
linestyles = ["-","--","-.",":"];
figure()
t = tiledlayout(1,5,'TileSpacing','Compact','Padding','Compact');
nexttile
plot(sigmas,E_Bws,'-*','Color',colors(1));
xlabel("$\sigma_{w}/a_{\rm{min}}$",'Interpreter','latex')
ylabel("mean($\left| {{{\hat B}_w}\left( \omega  \right) - \hat B\left( \omega  \right)} \right|$)",'Interpreter','latex')
title('(a)');
nexttile
plot(ratios,E_Bws2,'-*','Color',colors(1))
xlabel("$a_{\rm{max}}/a_{\rm{min}}$",'Interpreter','latex')
ylabel("mean($\left| {{{\hat B}_w}\left( \omega  \right) - \hat B\left( \omega  \right)} \right|$)",'Interpreter','latex')
title('(b)');
nexttile
plot(n_ps,E_Bws_n_p,'-*','Color',colors(1))
xlim([min(n_ps) max(n_ps)])
xlabel("$G$",'Interpreter','latex')
ylabel("mean($\left| {{{\hat B}_w}\left( \omega  \right) - \hat B\left( \omega  \right)} \right|$)",'Interpreter','latex')
title('(c)');
nexttile
plot(Ms,E_Bws_M,'-*','Color',colors(1))
xlabel("$M$",'Interpreter','latex')
ylabel("mean($\left| {{{\hat B}_w}\left( \omega  \right) - \hat B\left( \omega  \right)} \right|$)",'Interpreter','latex')
title('(d)');
load E_MQ;
nexttile
plot(Ms,E_Bws_MQ,'-*','Color',colors(1))
xticks(Ms(1:3:end))
ticks = strings(1,length(Ms(1:3:end)));
for i = 1:3
    ticks(i) = [num2str(Ms(i*3-2)) '(' num2str(Ns(i*3-2)) ')'];
end
xticklabels(ticks)
xlabel("$M(Q)$",'Interpreter','latex')
ylabel("mean($\left| {{{\hat B}_w}\left( \omega  \right) - \hat B\left( \omega  \right)} \right|$)",'Interpreter','latex')
title('(e)');

    