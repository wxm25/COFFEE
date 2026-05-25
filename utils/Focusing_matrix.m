function T = Focusing_matrix(theta, f, Omega, M, K, iter)
% Robust auto-focusing wideband DOA estimation, 2006 SP, Theorem IV-B and IV-C  
% theta is a row vector
N = length(Omega);
p = 2;
if iter == 1
    Spatial_Freq_interval = [-1/2,1/2];
    Q = focusing_integral(f*2, 1, Spatial_Freq_interval, Omega);
    for m = 1 : M
        [U,~,V] = svd(Q{m});
        T{m} = V*U';
    end
else
    Spatial_Freq_interval = sin(theta.')/2 + [-1/2/iter^p, 1/2/iter^p];  % modify
    Q = focusing_integral(f*2, 1, Spatial_Freq_interval, Omega);
    A0 = exp(1i*pi*kron(Omega-1,sin(theta)));
    for m = 1 : M
        A{m} = exp(1i*2*pi*f(m)*kron(Omega-1,sin(theta)));
        [U,~,V] = svd(A{m}*A0');
        U1 = U(:,1:K); V1 = V(:,1:K);
        X = orth(U(:,K+1:end)); Y = orth(V(:,K+1:end));
        [F,~,G] = svd(Y'*Q{m}'*X);
        T{m} = V1*U1' + Y*F*G'*X';
    end
end

end