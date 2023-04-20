%20230420 
%Generate Prob_G.mat to plot Fig 9 in this article. 
%@article{wang2023block, 
% title={Block-MUSIC in Blade Tip Timing: Performance Study of Block Snapshot Matrix}, 
% author={Wang, Zengkun and Yang, Zhibo and Wu, Shuming and Tian, Shaohua and Chen, Xuefeng}, 
% journal={Mechanical Systems and Signal Processing}, year={2023}}

clear all
close all

%parameter setting
omega = 6000/60; 
N_v = 24;             
A = [1 1];             
phase = [0 0];    
n_rev = 1000;     
sigma = 1;
N_mc = 200;
f = zeros(1,2);
f(1) = 873;       
delta_fs = [0.01:0.2:5];
normalized_f = delta_fs*1/N_v*0.01; 

%for different number of probes
n_ps = [1 2 4 8 16];
m = 80;
n_s = [129 65 33 17 9]; %number of snapshots
Prob_G = zeros(length(n_ps),length(delta_fs));
Rayleighs = 1./(m+(n_s-1).*n_ps);
for i_np = 1:length(n_ps)
    n_p = n_ps(i_np);
    waitbar(i_np/length(n_ps))
    n = n_s(i_np);
    for delta_f_i = 1:length(delta_fs)
        f(2) = 873+delta_fs(delta_f_i);
        delta = [0:15:(n_p-1)*15];
        %sampling time
        t = zeros(1,n_p*n_rev);
        for i = 0:n_rev-1
            for j = 1:length(delta)
                t(n_p*i+j)=(1/omega)*(i+delta(j)/360);
            end
        end
        %signal
        x = zeros(1,length(t));
        for i = 1:length(f)
            temp = A(i)*exp(1j*2*pi*f(i)*t+phase(i));
            x=x+temp;
        end
        temp_Es = zeros(1,N_mc);
        parfor n_mc = 1:N_mc
            % add noise    
            rng(2022+n_mc*i_np*delta_f_i);
            x_noi = x+sigma*randn(size(x));
            % 
            S2 = zeros(m,n);
            for i = 1:n
                S2(:,i)=x_noi((i-1)*n_p+1:(i-1)*n_p+m)';
            end
            [U2,D2,V2]=svd(S2);
            noise_space2 = U2(:,length(f)+1:end);

            w = [f mean(f)]; %
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
        Prob_G(i_np,delta_f_i) = mean(temp_Es);
    end
end

save("Prob_G.mat", "n_ps","normalized_f","Prob_G");

figure()
hold on
for i = 1:length(n_ps)
    plot(normalized_f,Prob_G(i,:));
end

    
    