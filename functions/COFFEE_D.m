function [DOA_rec, t_rec] = COFFEE_D(Y, N, Omega, M, L, K, f)
% Reference: X. Wu, Z. Yang, Z. Wei, R. Schober, and Z. Xu, "COFFEE: Covariance Fitting and Focusing for Wideband Direction-of-Arrival Estimation", IEEE Trans. Signal Processing, 2024.
% Written by Xunmeng Wu, 2024

I = eye(N); Gamma = I(Omega,:);
% Defaulat MM parameter setting
maxiter = 10; criterion = 10^(-3);
DOA_rec{1} = {};

for Iter = 1 : maxiter
    %     fprintf('Outer loop: %4g\n' ,Iter);
    Transform = Focusing_matrix(DOA_rec{Iter}, f, Omega, M, K, Iter);
    if L > N
        parfor m = 1 : M
            Y_Focus{m} = Transform{m} * Y{m};
            Rhat{m} = Y_Focus{m} * Y_Focus{m}' / L;
            [eig_vec, eig_val] = eig(Rhat{m});
            Yhat{m} = eig_vec * sqrtm(eig_val); 
            R_inv{m} = Gamma.'*pinv(Rhat{m})*Gamma;

            t_tmp = CVX(R_inv{m}, Yhat{m}, Omega, N);
            t(:,m) = t_tmp;
        end
    elseif L <= N  
        parfor m = 1 : M
            Y_Focus{m} = Transform{m} * Y{m};
            Rhat{m} = Y_Focus{m} * Y_Focus{m}' / L;
            [eig_vec, eig_val] = eigs(Rhat{m}*Rhat{m}, L, 'largestreal');  
            Yhat{m} = eig_vec * sqrtm(eig_val); 

            t_tmp = CVX_small_L(Yhat{m}, Omega, N, L);
            t(:,m) = t_tmp;
        end
    end

    t_rec{Iter+1} = t;
    [DOA, ~] = rootmusic(Toep(sum(t,2)),K,'corr');
    DOA = sort(asin(DOA/pi))'; DOA_rec{Iter+1} = DOA;

    % stopping criterion
    if Iter > 1
        res = norm(DOA_rec{Iter+1}-DOA_rec{Iter},2);
        if res < criterion
            break;
        end
    end
end

end

function t = CVX(R_inv, Yhat, Omega, N)
cvx_quiet true
cvx_solver sdpt3
cvx_begin sdp
variable t(2*N-1,1) complex;
variable Tt(N,N) hermitian;
variable W(N,N) hermitian;
variable Z(N,N) complex;
minimize abs(trace(R_inv*Tt + W)) + 2*sqrt(abs(trace(R_inv)))*norm(Yhat-Z(Omega,Omega),'fro');
subject to
Tt == Toep(t);
[W Z';
    Z Tt] == semidefinite(2*N,2*N);
cvx_end
end

function t = CVX_small_L(Yhat, Omega, N, L)
cvx_quiet true
cvx_solver sdpt3
cvx_begin sdp
variable t(2*N-1,1) complex;
variable Tt(N,N) hermitian;
variable W(L,L) hermitian;
variable Z(N,L) complex;
minimize abs(trace(Tt)) + abs(trace(W)) + 2*sqrt(N)*norm(Yhat-Z(Omega,:),'fro');
subject to
Tt == Toep(t);
[W Z';
    Z Tt] == semidefinite(N+L,N+L);
cvx_end
end

