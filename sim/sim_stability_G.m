%20230420 
%Generate E_G.mat plot Fig 7 in this article. 
%@article{wang2023block, 
% title={Block-MUSIC in Blade Tip Timing: Performance Study of Block Snapshot Matrix}, 
% author={Wang, Zengkun and Yang, Zhibo and Wu, Shuming and Tian, Shaohua and Chen, Xuefeng}, 
% journal={Mechanical Systems and Signal Processing}, year={2023}}

clear all
close all

%parameter setting
omega = 6000/60; %rotating frequency
N_v = 24;%number of virtual probes 
d_t = 1/omega/N_v;
f = [152 873];        %frequency 
A = [2 1];             %amplitude
phase = [0 0];    %phase
n_rev = 200;         %number of revolutions
q_2 = min(f(2:end)-f(1:end-1))*d_t;
q_1 = q_2*N_v;
m = 80; %length of snapshot: should be interger numbers of n_p
n = 80; %number of snapshot
sigma = 1;
n_ps = [4 8 10 16 20];
n_p = min(n_ps);
if q_2<1/n_p
    error(['q_2 = ' num2str(q_2) '<1/n_p, error!']);
end
up_lim_E_w = sigma*sqrt(2*m*log(2*m));
low_lim_sigma_K = min(A)*sqrt((m/n_p-1/q_1)*(n_p-1/q_2)*(n-1/q_1));
if up_lim_E_w>low_lim_sigma_K %verify that the norm of E is smaller than the smallest sigular value
    error(['sigma is too large']);
end

E_Bws_n_p = zeros(size(n_ps));
N_mc = 200;

for i_n_p = 1:length(n_ps)
    waitbar(i_n_p/length(n_ps))
    n_p = n_ps(i_n_p);
    temp_E = 0;
    delta = [0:15:(n_p-1)*15]; %probe position 
    %sampling time
    t = zeros(1,n_p*n_rev);
    for i = 0:n_rev-1
        for j = 1:length(delta)
            t(n_p*i+j)=(1/omega)*(i+delta(j)/360);
        end
    end

    parfor n_mc = 1:N_mc
        %signal generation
        x = zeros(1,length(t));
        for i = 1:length(f)
            temp = A(i)*exp(1j*2*pi*f(i)*t+phase(i));
            x=x+temp;
        end
        
        % noise-free case
        % snapshot matrix
        S = zeros(m,n);
        for i = 1:n
            S(:,i)=x((i-1)*n_p+1:(i-1)*n_p+m)';
        end
        [U1,D1,V1]=svd(S);
        noise_space = U1(:,length(f)+1:end);
        %steering vector
        w = 0.1:0.1:1200;
        noise_spectrum = zeros(size(w));
        steer_vectors = zeros(m,length(w));
        for j = 1:length(w)%
            steer_vectors(:,j) = exp(1j*t(1:m)*2*pi*w(j))'/sqrt(m);%normlized by sqrt(m)
        end
        
        % (music) traverse
        for i = 1:length(w)
            steer_vector = steer_vectors(:,i);
            noise_spectrum(i) = (((steer_vector'*noise_space)*(steer_vector'*noise_space)'));
        end
    
        % add noise    
        rng(2023+n_mc*i_n_p);
        x_noi = x+sigma*randn(size(x));
        % snapshot matrix
        S2 = zeros(m,n);
        for i = 1:n
            S2(:,i)=x_noi((i-1)*n_p+1:(i-1)*n_p+m)';
        end
        [U2,D2,V2]=svd(S2);
        noise_space2 = U2(:,length(f)+1:end);
        E = S2-S;
        % (music) traverse
        noise_spectrum2 = zeros(size(w));
        for i = 1:length(w)
            steer_vector = steer_vectors(:,i);
            noise_spectrum2(i) = (((steer_vector'*noise_space2)*(steer_vector'*noise_space2)'));
        end
        temp_E = temp_E+mean(abs(noise_spectrum-noise_spectrum2));
    end
    E_Bws_n_p(i_n_p) = temp_E/N_mc;
end

figure()
plot(n_ps,E_Bws_n_p);
save("E_G.mat", "n_ps","E_Bws_n_p");


    
    