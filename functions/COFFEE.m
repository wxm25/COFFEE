function [DOA_rec, t_rec, Iter_ADMM] = COFFEE(Y, N, Omega, M, L, K, f)
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
        for m = 1 : M
            Y_Focus{m} = Transform{m} * Y{m};
            Rhat{m} = Y_Focus{m} * Y_Focus{m}' / L;
            [eig_vec, eig_val] = eig(Rhat{m});
            Yhat{m} = eig_vec * sqrtm(eig_val);  
            R_inv{m} = Gamma.'*pinv(Rhat{m})*Gamma;
        end
        % Default ADMM parameter setting
        param.tol_abs = 1e-2;   % absolute tolerance
        param.tol_rel = 1e-3;   % relevant tolerance
        param.verbose = 0; % set to 1 for printing, 0 for no output
        param.mu1 = ones(M,1)*1;          % penalty parameter
        param.mu2 = 1;          % penalty parameter
        param.maxiter = 1000;   % maximum iterations 
        param.alpha = 2;
        param.gamma = 10;        % tune the penalty parameter
        if Iter == 2 
            param.tol_abs = 1e-3;   
            param.tol_rel = 1e-4;
        elseif Iter > 2 
            param.tol_abs = 1e-4;
            param.tol_rel = 1e-5;
        end
        [t, iter] = ADMM(R_inv, Yhat, N, M, K, param);
    elseif L <= N 
        for m = 1 : M
            Y_Focus{m} = Transform{m} * Y{m};
            Rhat{m} = Y_Focus{m} * Y_Focus{m}' / L;
            [eig_vec, eig_val] = eigs(Rhat{m}*Rhat{m}, L, 'largestreal');  %%
            Yhat{m} = eig_vec * sqrtm(eig_val);  
        end
        % Default ADMM parameter setting
        param.tol_abs = 1e-2;   % absolute tolerance
        param.tol_rel = 1e-3;   % relevant tolerance
        param.verbose = 0; % set to 1 for printing, 0 for no output
        param.mu1 = ones(M,1)*1;          % penalty parameter
        param.mu2 = 1;          % penalty parameter
        param.maxiter = 1000;   % maximum iterations 
        param.alpha = 2;
        param.gamma = 10;        % tune the penalty parameter
        if Iter == 2 
            param.tol_abs = 1e-3;   % absolute tolerance
            param.tol_rel = 1e-4;
        elseif Iter > 2 
            param.tol_abs = 1e-4;   % absolute tolerance
            param.tol_rel = 1e-5;
        end
        [t, iter] = ADMM_small_L(Rhat, Yhat, N, M, K, L, param);
    end

    Iter_ADMM(Iter) = iter;
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
