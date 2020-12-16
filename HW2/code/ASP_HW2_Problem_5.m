clear, clc
mkdir figure

load('ASP_Problem_5.mat')
x = filter([1 1], [1 -1/6 -1/6], v);
[~, n] = size(x);

k = [0:100].';
r = (108 / 35) * (2) .^ (-abs(-k)) - (18 / 35) * (-3) .^ (-abs(k));
[a, P, kappa] = ASP_Levinson_Durbin(r(1:2));

% Problem 5 (b)
fM_1 = a{2}'*[x; rshift(x, 1)];

figure;
plot(1:n, real(fM_1), 'r-');
hold on;
plot(1:n, imag(fM_1), 'b-');
xlabel('Time index');
ylabel('Amplitude');
title('Real and Imag part of f_1(n)');
legend({'Real part', 'Imag part'})
grid on;
saveas(gcf, 'figure/ASP_HW2_Problem_5b.jpg');

Pf_1 = fM_1 * fM_1' / n;
disp('The average power of f1(n) is:');
disp(Pf_1);

% Problem 5 (c)
bM_1 = flip(a{2})'*[rshift(x, -1); x];
x_abs = abs(x);
xr_abs = abs(rshift(x, -1));

figure;
plot(1:n, real(bM_1), 'r-');
hold on;
plot(1:n, imag(bM_1), 'b-');
xlabel('Time index');
ylabel('Amplitude');
title('Real and Imag part of b_1(n)');
legend({'Real part', 'Imag part'})
grid on;
saveas(gcf, 'figure/ASP_HW2_Problem_5c.jpg');

Pb_1 = bM_1 * bM_1' / n;
disp('The average power of b1(n) is:');
disp(Pb_1);

% Problem 5 d
n_iter = 10;
Pf = [];
Pb = [];
x_forwards = x;
% x_backwards = x;
x_backwards = rshift(x, n_iter);
for m = 1 : n_iter
    x_forwards = [x_forwards; rshift(x, m)];
    % x_backwards = [rshift(x, -m); x_backwards];
    x_backwards = [rshift(x, n_iter-m); x_backwards];
end
for m = 1 : n_iter
    [a, P, kappa] = ASP_Levinson_Durbin(r(1:(m+1)));
    % forward error
    
    fM = a{m+1}'*x_forwards(1:m+1, :);
    Pf = [Pf fM*fM'/n];
    % backward error
    % bM = flip(a{m+1})'*x_backwards(end-m:end, :);
    bM = flip(a{m+1})'*x_backwards(1:m+1, :);
    Pb = [Pb bM*bM'/n];
end

H = @(z) ( (1+z.^(-1)) ./ (1-1/2*z.^(-1)) ./ (1+1/3*z.^(-1)) );
S_x = @(z) ( H(z) .* conj(H(1./conj(z))) );
func = @(x) (log(S_x(exp(2*1i*pi*x))));
lower_bound = exp(integral(func, -1/2, 1/2));

figure;
plot(1:n_iter, Pf, 'r-^', 'MarkerSize', 5,'LineWidth', 2);
hold on;
plot(1:n_iter, Pb ,'b-o', 'MarkerSize', 5,'LineWidth', 2);
plot(1:n_iter, P(2:end), 'g-o', 'MarkerSize', 5,'LineWidth', 2);
plot(1:n_iter, lower_bound*ones(1,n_iter), 'k-o', 'MarkerSize', 5,'LineWidth', 2);
xlabel('Index');
ylabel('Error power');
title('Error power - Index');
legend({'Forward Prediction Error', 'Backward Prediction Error', 'Prediction Error', 'Lower bound'});
grid on;
saveas(gcf, 'figure/ASP_HW2_Problem_5d.jpg');


function x_out = rshift(x_in, n_shift)
    if n_shift >= 0
        x_out = circshift(x_in, n_shift);
        x_out(1:n_shift) = 0;
    else
        x_out = circshift(x_in, n_shift);
        x_out(end+n_shift+1:end) = 0;
    end
end
