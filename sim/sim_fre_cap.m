%20230420 frequency capacity
%Fig 6 in this article. 
%@article{wang2023block, 
% title={Block-MUSIC in Blade Tip Timing: Performance Study of Block Snapshot Matrix}, 
% author={Wang, Zengkun and Yang, Zhibo and Wu, Shuming and Tian, Shaohua and Chen, Xuefeng}, 
% journal={Mechanical Systems and Signal Processing}, year={2023}}

clear all
close all

rng(2023)
%parameter setting
omega = 6000/60; %rotating speed
f = sort(round(1200*unifrnd(0,1,[1,7])));%frequency
A = [1 1 1 1 1 1 1 1 1 1];             %amplitude
phase = [0 0 0 0 0 0 0 0 0 0];    
n_rev = 1000 ;         %number of revolutions

%% n_p = 1
delta = [0]; %propbe position
n_p = length(delta);   
%signal generation
t = zeros(1,n_p*n_rev);
for i = 0:n_rev-1
    for j = 1:length(delta)
        t(n_p*i+j)=(1/omega)*(i+delta(j)/360);
    end
end
x = zeros(1,length(t));
for i = 1:length(f)
    temp = A(i)*exp(1j*2*pi*f(i)*t+phase(i));
    x=x+temp;
end

%snapshot matrix1
m = 11*n_p; 
n = 11; 
S = zeros(m,n);
for i = 1:n
    S(:,i)=x((i-1)*n_p+1:(i-1)*n_p+m)';
end
S=S/sqrt(m*n);
[U1,D1,V1]=svd(S);
sin_values_12_1np = diag(D1);

%snapshot matrix2
m2 = 7*n_p; 
n2 = 7; 
S2 = zeros(m2,n2);
for i = 1:n2
    S2(:,i)=x((i-1)*n_p+1:(i-1)*n_p+m2)';
end
[U2,D2,V2]=svd(S2);
sin_values_10_1np = diag(D2);


%% n_p = 2
delta = [0 15]; 
n_p = length(delta);   
%signal generation
t = zeros(1,n_p*n_rev);
for i = 0:n_rev-1
    for j = 1:length(delta)
        t(n_p*i+j)=(1/omega)*(i+delta(j)/360);
    end
end
x = zeros(1,length(t));
for i = 1:length(f)
    temp = A(i)*exp(1j*2*pi*f(i)*t+phase(i));
    x=x+temp;
end

%snapshot matrix1
m = 11*n_p; 
n = 11; 
S = zeros(m,n);
for i = 1:n
    S(:,i)=x((i-1)*n_p+1:(i-1)*n_p+m)';
end
[U1,D1,V1]=svd(S);
sin_values_12_2np = diag(D1);

%snapshot matrix2
m2 = 7*n_p; 
n2 = 7; 
S2 = zeros(m2,n2);
for i = 1:n2
    S2(:,i)=x((i-1)*n_p+1:(i-1)*n_p+m2)';
end
S2 = S2/sqrt(m2*n2);
[U2,D2,V2]=svd(S2);
sin_values_10_2np = diag(D2);

% plot1
colors = ["#1A2B3E","#185F8C","#97BEA1","#A4A750","#8B5C0C"];
linestyles = ["-","--","-.",":"];
figure()
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
nexttile

semilogy([1:length(sin_values_10_1np)],sin_values_10_1np,...
    linestyles(1),'Color',colors(1));
hold on
title('(a)');
semilogy([1:length(sin_values_12_1np)],sin_values_12_1np, ...
    linestyles(2),'Color',colors(2));
semilogy([1:length(sin_values_10_2np)],sin_values_10_2np, ...
    linestyles(3),'Color',colors(3));
semilogy([1:length(sin_values_12_2np)],sin_values_12_2np, ...
    linestyles(4),'Color',colors(4));
xlim([1,11]);
ylim([0.99*min(sin_values_12_1np),1.01*max(sin_values_12_2np)]);
xlabel("Index")
ylabel("Singular value")
hl = legend("$G=1,L=7$","$G=1,L=11$","$G=2,L=7$",...
    "$G=2,L=11$",'Interpreter','latex','location','south','numcolumns',1);
set(hl,'Box','off');
box on


%% the figure versus noise 
delta = [0 15]; 
n_p = length(delta);   
%signal generation
t = zeros(1,n_p*n_rev);
for i = 0:n_rev-1
    for j = 1:length(delta)
        t(n_p*i+j)=(1/omega)*(i+delta(j)/360);
    end
end
x = zeros(1,length(t));
for i = 1:length(f)
    temp = A(i)*exp(1j*2*pi*f(i)*t+phase(i));
    x=x+temp;
end

sigma_ws = [0.2:0.2:0.8];
m = 50*n_p; 
n = 50; 
sin_values = zeros(n,length(sigma_ws));
for i_sigma_w = 1:length(sigma_ws)
    sigma = sigma_ws(i_sigma_w);
    noise = sigma*randn(size(x));
    var = sum(noise*noise')/length(noise);
    x_temp = x+noise;
    %snapshot matrix
    S = zeros(m,n);
    for i = 1:n
        S(:,i)=x_temp((i-1)*n_p+1:(i-1)*n_p+m)';
    end
    S = S/sqrt(n*m);
    [U1,D1,V1]=svd(S);
    sin_values(:,i_sigma_w) = diag(D1);
end


% plot2
colors = ["#1A2B3E","#185F8C","#97BEA1","#A4A750","#8B5C0C"];
linestyles = ["-","--","-.",":"];
nexttile
title('(b)');
for i = 1:length(sigma_ws)
    hold on
    plot([1:length(sin_values(:,i))],sin_values(:,i),...
        linestyles(i),'Color',colors(i));
end
xlim([1,25]);
ylim([0.99*min(min(sin_values)),1.01*max(max(sin_values))]);
xlabel("Index")
ylabel("Singular value")
hl = legend("$\sigma_w=0.2$","$\sigma_w=0.4$","$\sigma_w=0.6$",...
    "$\sigma_w=0.8$",'Interpreter','latex','location','south','numcolumns',1);
set(hl,'Box','off');
box on

    
    