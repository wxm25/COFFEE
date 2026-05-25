function Q = Proj_using_SVD(Input, K)
if nargin < 2
    [p, q] = eig((Input + Input') / 2);
    dq = real(diag(q));
    idxpos = (dq > 0);
    Q = p(:,idxpos) * diag(dq(idxpos)) * p(:,idxpos)';
elseif nargin == 2
    [U,D,V] = svd(Input);
    d = diag(D);
    Q = U(:,1:K) * diag(d(1:K)) * V(:,1:K)';
end

end

