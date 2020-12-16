M = 5;
FIR_len = 30;

h = @(n) ( (18/25)*(1/2).^n + (7/25)*(-1/3).^n );
h = h([0:FIR_len-1]);

R = zeros(M, M);
for ii = 1 : M
    for jj = 1 : M
        R(ii, jj) = h*rshift(h,jj-ii)';
    end
end
p = [1 ; zeros(M-1, 1)];
w = inv(R)*p;
J_min = 1 - p'*w;

function x_out = rshift(x_in, n_shift)
    if n_shift >= 0
        x_out = circshift(x_in, n_shift);
        x_out(1:n_shift) = 0;
    else
        x_out = circshift(x_in, n_shift);
        x_out(end+n_shift+1:end) = 0;
    end
end
