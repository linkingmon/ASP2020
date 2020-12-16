function J = ASP_NLMS(x, d, M, mu_NLMS)
    % x = [x(1) ; x(2) ; ... ]
    % d = [d(1) ; d(2) ; ... ]
    N = size(x, 1);
    w = zeros(M, 1);
    x_extend = [zeros(M-1, 1); x];
    e = zeros(N, 1);
    y = zeros(N, 1);
    J = zeros(N, 1);

    for idx = 1:N
        x = flip(x_extend(idx:idx+M-1));
        y(idx, 1) = w' * x;
        e(idx, 1) = d(idx, 1) - y(idx, 1);
        J(idx, 1) = e(idx, 1)' * e(idx, 1);
        mu = mu_NLMS / (x'*x);
        w = w + mu * x * conj(e(idx));
    end
end