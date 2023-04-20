% 20230420 
% A simple illustration of Block-MUSIC 
% Please cite
%@article{wang2023block, 
% title={Block-MUSIC in Blade Tip Timing: Performance Study of Block Snapshot Matrix}, 
% author={Wang, Zengkun and Yang, Zhibo and Wu, Shuming and Tian, Shaohua and Chen, Xuefeng}, 
% journal={Mechanical Systems and Signal Processing}, year={2023}}
%@article{wang2020improved, 
% title={An improved multiple signal classification for nonuniform sampling in blade tip timing}, 
% author={Wang, Zengkun and Yang, Zhibo and Wu, Shuming and Li, Haoqi and Tian, Shaohua and Chen, Xuefeng}, 
% journal={IEEE Transactions on Instrumentation and Measurement}, 
% year={2020}}

clear all
close all
rng(2023);

%parameter setting
omega = 6000/60; %rotating frequency
f = [152 873];        %frequency
A = [1.2 1];          %amplitude
phase = [0 0];        %phase
n_rev = 500;          %number of revolutions
sigma = 0.5;          %standard variance of noise

delta = [0 15.5 32.1 47.2]; %install angles of probes 
n_p = length(delta); 

%signal generation
t = zeros(1,n_p*n_rev); %sampling time
for i = 0:n_rev-1
    for j = 1:length(delta)
        t(n_p*i+j)=(1/omega)*(i+delta(j)/360);
    end
end
x = zeros(1,length(t)); %vibration displacement
for i = 1:length(f)
    temp = A(i)*sin(2*pi*f(i)*t+phase(i));
    x=x+temp;
end

% add_noise    
x_noi = x+sigma*randn(size(x));

%Block-MUSIC
m = 30*n_p; %snapshot length
n = m;      %snapshot number
S = zeros(m,n); 
for i = 1:n
    S(:,i)=x_noi((i-1)*n_p+1:(i-1)*n_p+m)';
end
[U,D,V]=svd(S);
evs = diag(D);
K_max = 10;% max possible number of frequencies
for i = 1:length(evs)
    if sum(evs(1:i))/sum(evs)>0.95
        noise_space = U(:,i+1:end);
        break
    end
    if i==length(evs)
        noise_space = U(:,2*K_max+1:end);
    end
end

w = 0.1:0.1:1200;
steer_vectors = zeros(m,length(w));
for j = 1:length(w)
    steer_vectors(:,j) = exp(1j*t(1:m)*2*pi*w(j))'/sqrt(m);
end

% Frequency traverse
noise_spectrum = zeros(size(w));
for i = 1:length(w)
    steer_vector = steer_vectors(:,i);
    noise_spectrum(i) = 1./(((steer_vector'*noise_space)*(steer_vector'*noise_space)'));
end

figure()
plot(w,noise_spectrum);
xlabel("frequency/Hz")
ylabel("pseudo amplitude")

    
    