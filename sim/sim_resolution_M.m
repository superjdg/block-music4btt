%20230420 
%plot Fig 8 in this article. 
%@article{wang2023block, 
% title={Block-MUSIC in Blade Tip Timing: Performance Study of Block Snapshot Matrix}, 
% author={Wang, Zengkun and Yang, Zhibo and Wu, Shuming and Tian, Shaohua and Chen, Xuefeng}, 
% journal={Mechanical Systems and Signal Processing}, year={2023}}

clear all
close all

%parameter setting
omega = 6000/60; %rotating frequency
delta = [0 15 30 45]; %probe position
n_p = length(delta);  %probe number
N_v = 24;             %virtual probe number 
A = [1 1];             %amplitude
phase = [0 0];   
n_rev = 1000;         
sigma = 1;
N_mc = 200;
f = zeros(1,2);
f(1) = 873;
delta_fs = [0.01:0.2:5];
normalized_f = delta_fs*1/N_v*0.01; 

%for different snapshot lengths
Ms = [64:4:76];
n = 34; %number of snapshots
Rayleigh = 1./(Ms+(n-1)*n_p);
P_rec_c = zeros(length(Ms),length(delta_fs));
for delta_f_i = 1:length(delta_fs)
    waitbar(delta_f_i/length(delta_fs))
    f(2) = 873+delta_fs(delta_f_i);
    %sampling time
    t = zeros(1,n_p*n_rev);
    for i = 0:n_rev-1
        for j = 1:length(delta)
            t(n_p*i+j)=(1/omega)*(i+delta(j)/360);
        end
    end
    %signal generation
    x = zeros(1,length(t));
    for i = 1:length(f)
        temp = A(i)*exp(1j*2*pi*f(i)*t+phase(i));
        x=x+temp;
    end

    for i_M = 1:length(Ms)
        
        m = Ms(i_M);
        temp_Es = zeros(1,N_mc);
        parfor n_mc = 1:N_mc
            % add_noise    
            rng(2023+n_mc*i_M*delta_f_i);
            x_noi = x+sigma*randn(size(x));
            % snapshot
            S2 = zeros(m,n);
            for i = 1:n
                S2(:,i)=x_noi((i-1)*n_p+1:(i-1)*n_p+m)';
            end
            [U2,D2,V2]=svd(S2);
            noise_space2 = U2(:,length(f)+1:end);

            w = [f mean(f)]; 
            noise_spectrum = zeros(size(w));
            steer_vectors = zeros(m,length(w));
            for j = 1:length(w)
                steer_vectors(:,j) = exp(1j*t(1:m)*2*pi*w(j))'/sqrt(m);%normlized by sqrt(m)
            end
            % (music)
            noise_spectrum2 = zeros(size(w));
            for i = 1:length(w)
                steer_vector = steer_vectors(:,i);
                noise_spectrum2(i) = (((steer_vector'*noise_space2)*(steer_vector'*noise_space2)'));
            end
            if noise_spectrum2(3)>=noise_spectrum2(1)... 
                && noise_spectrum2(3)>=noise_spectrum2(2)
                temp_Es(n_mc) = 1;
            end
        end
        P_rec_c(i_M,delta_f_i) = mean(temp_Es);
    end
end



%for different snapshot numbers
ns = [34:4:46];
m = 64; %length of snapshot
Rayleigh_n = 1./(m+(ns-1)*n_p);
P_rec_n = zeros(length(ns),length(delta_fs));
for delta_f_i = 1:length(delta_fs)
    waitbar(delta_f_i/length(delta_fs))
    f(2) = 873+delta_fs(delta_f_i);
    %sampling time 
    t = zeros(1,n_p*n_rev);
    for i = 0:n_rev-1
        for j = 1:length(delta)
            t(n_p*i+j)=(1/omega)*(i+delta(j)/360);
        end
    end
    %signal generation
    x = zeros(1,length(t));
    for i = 1:length(f)
        temp = A(i)*exp(1j*2*pi*f(i)*t+phase(i));
        x=x+temp;
    end

    for i_n = 1:length(ns)
        
        n = ns(i_n);
        temp_Es = zeros(1,N_mc);
        parfor n_mc = 1:N_mc
            % add noise    
            rng(2023+n_mc*i_n*delta_f_i*10);
            x_noi = x+sigma*randn(size(x));
            % snapshot matrix
            S3 = zeros(m,n);
            for i = 1:n
                S3(:,i)=x_noi((i-1)*n_p+1:(i-1)*n_p+m)';
            end
            [U3,D3,V3]=svd(S3);
            noise_space3 = U3(:,length(f)+1:end);

            w = [f mean(f)]; 
            noise_spectrum = zeros(size(w));
            steer_vectors = zeros(m,length(w));
            for j = 1:length(w)
                steer_vectors(:,j) = exp(1j*t(1:m)*2*pi*w(j))'/sqrt(m);%normlized by sqrt(m)
            end
            % (music) 
            noise_spectrum3 = zeros(size(w));
            for i = 1:length(w)
                steer_vector = steer_vectors(:,i);
                noise_spectrum3(i) = (((steer_vector'*noise_space3)*(steer_vector'*noise_space3)'));
            end
            if noise_spectrum3(3)>=noise_spectrum3(1)... 
                && noise_spectrum3(3)>=noise_spectrum3(2)
                temp_Es(n_mc) = 1;
            end
        end
        P_rec_n(i_n,delta_f_i) = mean(temp_Es);
    end
end

%plot
colors = ["#1A2B3E","#185F8C","#97BEA1","#A4A750","#8B5C0C"];
linestyles = ["-","--","-.",":"];
t = tiledlayout(1,2,'TileSpacing','Compact','Padding','Compact');
nexttile
hold on
for i = 1:size(P_rec_c,1)
    plot(normalized_f,P_rec_c(i,:),linestyles(i),'Color',colors(i));
end
xlim([min(normalized_f),max(normalized_f)]);
ylim([0, 1]);
box on
hl = legend("$M=64$","$M=68$","$M=72$",...
    "$M=76$",'Interpreter','latex','location','north','numcolumns',1);
set(hl,'FontName','Times New Roman','FontSize',6,'FontWeight','normal')
set(hl,'Box','off');
xlabel("$\omega_r-\omega_l$",'Interpreter','latex')
ylabel("Probability",'Interpreter','latex');
title('(a)');

nexttile
hold on
for i = 1:size(P_rec_n,1)
    plot(normalized_f,P_rec_n(i,:),linestyles(i),'Color',colors(i));
end
xlim([min(normalized_f),max(normalized_f)]);
ylim([0, 1]);
box on
hl = legend("$Q=34$","$Q=38$","$Q=42$",...
    "$Q=46$",'Interpreter','latex','location','north','numcolumns',1);
set(hl,'FontName','Times New Roman','FontSize',6,'FontWeight','normal')
set(hl,'Box','off');
xlabel("$\omega_r-\omega_l$",'Interpreter','latex')
ylabel("Probability");
title('(b)');
set(gcf,'unit','centimeters','position',[3 5 8 6])

    
    