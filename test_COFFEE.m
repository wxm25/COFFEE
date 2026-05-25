%% Wideband DOA Estimation via COFFEE %%
warning off
addpath('utils'); addpath('functions'); addpath('D:\Matlab2019b\cvx')
clear all; clc; close all;

%% Signal generation
N = 11; Omega = [1:N]'; L = 100; 
M0 = 64; M = 11;
delta = 0.2;
K = 2; w_theta = [-0.1, -0.1+delta/N]; theta = asin(2*w_theta);
Index = [(M0-M+3)/2:(M0+M+1)/2];
Total_signal_power = 0;
for m = 1 : M
    f(m) = (Index(m)-1)/M0;
    A{m} = exp(1i*2*pi*f(m)*kron(Omega-1,sin(theta)));
    cor = eye(K);
    am = 1 + abs(randn(1)); pow = am*am'.*diag(sqrt([1:K]));
    P{m} = pow.*cor;
    [U,D] = eig(P{m});
    L_matrix = U*sqrt(D);
    S{m} = L_matrix * (randn(K,L)+1i*randn(K,L))/sqrt(2);
    Total_signal_power = Total_signal_power + sum(sum(P{m}));
    X{m} = A{m} * S{m};
end
SNR = 20;
Total_sigma = Total_signal_power / 10^(SNR/10);
beta = rand(M,1);
sigma = beta / sum(beta) * Total_sigma;
for m = 1 : M
    noise{m} = sqrt(sigma(m)*eye(N))*(randn(N,L) + 1i*randn(N,L))/sqrt(2);
    Y{m} = X{m} + noise{m};
end

%% COFFEE
t1 = clock; [theta_C, t_C] = COFFEE(Y, N, Omega, M, L, K, f); t2 = clock;
MSE_C = norm(theta_C{end}-theta,2)^2; Time_C = etime(t2,t1);

%% a decoupled version of COFFEE
t1 = clock; [theta_D, t_D] = COFFEE_D(Y, N, Omega, M, L, K, f); t2 = clock;
MSE_D = norm(theta_D{end}-theta,2)^2; Time_D = etime(t2,t1);


