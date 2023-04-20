%20230420 
%Generate E_sigma_amin.mat plot Fig 7 in this article. 
%Refer to sim_stability_G.m for the variable meanings 
%@article{wang2023block, 
% title={Block-MUSIC in Blade Tip Timing: Performance Study of Block Snapshot Matrix}, 
% author={Wang, Zengkun and Yang, Zhibo and Wu, Shuming and Tian, Shaohua and Chen, Xuefeng}, 
% journal={Mechanical Systems and Signal Processing}, year={2023}}
clear all
close all

%parameter
omega = 6000/60; 
delta = [0 15 30 45]; 
n_p = length(delta);   
N_v = 24;
d_t = 1/omega/N_v;
f = [152 873];        
q_2 = min(f(2:end)-f(1:end-1))*d_t;
q_1 = q_2*N_v;
if q_2<1/n_p
    error(['q_1 = ' num2str(q_1) '<1/n_p, error!']);
end


A = [2 1];            
phase = [0 0];    
n_rev = 1000;        

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

%% noise-free
% snapshot matrix
m = 20*n_p; 
n = 80; 
S = zeros(m,n);
for i = 1:n
    S(:,i)=x((i-1)*n_p+1:(i-1)*n_p+m)';
end
[U1,D1,V1]=svd(S);
noise_space = U1(:,length(f)+1:end);

w = 0.1:0.1:1200;
noise_spectrum = zeros(size(w));
steer_vectors = zeros(m,length(w));
for j = 1:length(w)
    steer_vectors(:,j) = exp(1j*t(1:m)*2*pi*w(j))'/sqrt(m);%normlized by sqrt(m)
end

% (music) 
for i = 1:length(w)
    steer_vector = steer_vectors(:,i);
    noise_spectrum(i) = (((steer_vector'*noise_space)*(steer_vector'*noise_space)'));
end

%% noisy case
sigmas = [0:0.1:1];
E_Bws = zeros(size(sigmas));
N_mc = 200;

for i_sigma = 1:length(sigmas)
    waitbar(i_sigma/length(sigmas))
    sigma = sigmas(i_sigma);
    temp_E = 0;
    up_lim_E_w = sigma*sqrt(2*m*log(2*m));
    low_lim_sigma_K = min(A)*sqrt((m/n_p-1/q_1)*(n_p-1/q_2)*(n-1/q_1));
    if up_lim_E_w>low_lim_sigma_K 
        break
    end
    if sigma==0
       E_Bws(i_sigma) = temp_E;
       continue
    end
    parfor n_mc = 1:N_mc
        rng(2023+n_mc*i_sigma);
        x_noi = x+sigma*randn(size(x));
        % 
        S2 = zeros(m,n);
        for i = 1:n
            S2(:,i)=x_noi((i-1)*n_p+1:(i-1)*n_p+m)';
        end
        [U2,D2,V2]=svd(S2);
        noise_space2 = U2(:,length(f)+1:end);
        E = S2-S;
        % (music) 
        noise_spectrum2 = zeros(size(w));
        for i = 1:length(w)
            steer_vector = steer_vectors(:,i);
            noise_spectrum2(i) = (((steer_vector'*noise_space2)*(steer_vector'*noise_space2)'));
        end
        temp_E = temp_E+mean(abs(noise_spectrum-noise_spectrum2));
    end
    E_Bws(i_sigma) = temp_E/N_mc;
end

% figure()
% plot(sigmas,E_Bws);
save("E_sigma_amin.mat", "sigmas","E_Bws");

    
    