function [t_output, iter] = ADMM(R_inv, Yhat, N, M, K, param)
% Reference: X. Wu, Z. Yang, Z. Wei, R. Schober, and Z. Xu, "COFFEE: Covariance Fitting and Focusing for Wideband Direction-of-Arrival Estimation", IEEE Trans. Signal Processing, 2024.
% Written by Xunmeng Wu, 2024

T0 = clock;
% initialization
tol_abs = param.tol_abs;
tol_rel = param.tol_rel;
verbose = param.verbose;
mu1 = param.mu1;
mu2 = param.mu2;
maxiter = param.maxiter;
alpha = param.alpha;
gamma = param.gamma;
for m = 1 : M
    Z{m} = Yhat{m};
    Tt{m} = Yhat{m}*Yhat{m}';  
    W{m} = zeros(N,N);
    Lambda_m{m} = zeros(2*N,2*N);
end
Lambda = zeros(N,N);
Tt_sum = sum(cat(3,Tt{:}),3)/M;

if verbose
    fprintf('iter | res_prim   target    | res_dual  target    | iter     total\n')
    for k=1:68, fprintf('-'); end
    fprintf('\n');
end

for iter = 1 : maxiter
    T1 = clock;

    % update Q & Q{m}
    Q = Proj_using_SVD(Tt_sum - Lambda/mu2, K);  
    for m = 1 : M
        Big{m} = [W{m}, Z{m}'; Z{m}, Tt{m}];
        Q_m{m} = Proj_using_SVD(Big{m} - Lambda_m{m}/mu1(m));
    end

    % update W, Z, t
    for m = 1 : M
        W_new{m} = Q_m{m}(1:N,1:N) + Lambda_m{m}(1:N,1:N)/mu1(m) - eye(N)/mu1(m);
        L_matrix = Yhat{m} - Q_m{m}(N+1:2*N,1:N) - Lambda_m{m}(N+1:2*N,1:N)/mu1(m);
        if norm(L_matrix,'fro') < 1e-6
            Z_new{m} = Yhat{m};
        else
            Z_new{m} = Yhat{m} - max(1 - sqrt(trace(R_inv{m}))/(mu1(m)*norm(L_matrix,'fro')), 0) * L_matrix;
        end
        tmp{m} = R_inv{m}/mu1(m) - Q_m{m}(N+1:2*N,N+1:2*N) - Lambda_m{m}(N+1:2*N,N+1:2*N)/mu1(m);
    end

    tmp_sum = sum(cat(3,tmp{:}),3);
    for m = 1 : M
        t_new{m} = toeplitz_approx( (M*(Q+Lambda/mu2) + tmp_sum)/(mu1(m)/mu2*M+1)/M - tmp{m} );    
        Tt_new{m} = Toep(t_new{m});
        Big_new{m} = [W_new{m}, Z_new{m}'; Z_new{m}, Tt_new{m}];
    end

    % update Lambda & Lambda_m
    Tt_new_sum = sum(cat(3,Tt_new{:}),3)/M;
    Lambda = Lambda + mu2*(Q - Tt_new_sum);
    for m = 1 : M
        Lambda_m{m} = Lambda_m{m} + mu1(m)*(Q_m{m} - Big_new{m});
    end

    % stopping criterion (Boyd's rule)
    for m = 1 : M
        err_prim1(m) = norm(Big_new{m} - Q_m{m}, 'fro');
        Big_diff{m} = Big_new{m} - Big{m};
        err_dual1(m) = norm( mu1(m)*[vec(Big_diff{m}(1:N,1:N)); vec(Big_diff{m}(N+1:end,1:N)); toeplitz_approx(Big_diff{m}(N+1:end,N+1:end))],'fro' );
        tol_prim1(m) = sqrt((2*N)*(2*N))*tol_abs + tol_rel*max(norm(Big_new{m},'fro'), norm(Q_m{m},'fro')); 

        dual_var_adj(m) = norm([vec(Lambda_m{m}(1:N,1:N)); vec(Lambda_m{m}(N+1:end,1:N)); toeplitz_approx(Lambda_m{m}(N+1:end,N+1:end))],'fro');
        tol_dual1(m) = sqrt((N^2+N^2+2*N-1)*1)*tol_abs + tol_rel*dual_var_adj(m);  
        
        test(m) = err_prim1(m) < tol_prim1(m) && err_dual1(m) < tol_dual1(m);
    end

    % Q, Tt_sum
    err_prim2 = norm(Tt_new_sum - Q,'fro');
    Tt_sum_diff = Tt_new_sum - Tt_sum; err_dual2 = norm( mu2*toeplitz_approx(Tt_sum_diff),'fro');
    tol_prim2 = sqrt(N^2)*tol_abs + tol_rel*max(norm(Tt_new_sum,'fro'), norm(Q,'fro'));
    tol_dual2 = sqrt(N^2)*tol_abs + tol_rel*norm(vec(Lambda),'fro');

    if verbose
        fprintf('%4g | %.3d %.3d | %.3d %.3d | %.2d %.2d\n',...
            iter,err_prim2,tol_prim2,err_dual2,tol_dual2,etime(clock,T1),sum(test));  %etime(clock,T0)
    end
    if err_prim2 < tol_prim2 && err_dual2 < tol_dual2 && sum(test) == M
        if verbose,  fprintf('\n'); end
        break;
    end
    
    % a varying penalty parameter scheme 
    % fix mu1 since it corresponds to the convex constraint
    if err_prim2 > gamma * err_dual2
        mu2 = alpha * mu2;
    elseif err_prim2 < err_dual2 / gamma
        mu2 = mu2 / alpha;
    end

    W = W_new; Z = Z_new; Tt = Tt_new; Tt_sum = Tt_new_sum;
end

for m = 1 : M
    t_output(:,m) = t_new{m};
end

end