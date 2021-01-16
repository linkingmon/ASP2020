mkdir figure
disp('Loading numerical values the Kalman filter ...')
load('ASP_HW4_Problem_5.mat')
[N, L] = size(Y_tilde);
[~, M] = size(Q1_n);
disp(['M = ', int2str(M)]);
disp(['N = ', int2str(N)]);
disp(['L = ', int2str(L)]);

disp('Running Kalman filter ...')
x_tilde_i = zeros(M, 1);
K_i = eye(M);
x_tilde_o = zeros(M, L);
for ii = 1 : L
    y_tilde_i = Y_tilde(:, ii);
    [K_o, x_n1_tilde_o, x_tilde_o(:,ii)] = KalmanFilter(F_n1_n, C_n, Q1_n, Q2_n, K_i, x_tilde_i, y_tilde_i);
    K_i = K_o;
    x_tilde_i = x_n1_tilde_o;
end

disp('(i) Plotting real part of x ...')
figure;
for ii = 1 : M
    subplot(M,1,ii)
    plot(1:L, real(x_tilde_o(ii,:)))
    xlabel('Time index (n)')
    ylabel('Value')
    title(['Real part of x_', int2str(ii), '(n | Y_n)'])
    xlim([1 L])
    grid on
end
saveas(gcf, ['figure/ASP_HW4_Problem_5_Real', '.jpg']);
saveas(gcf, ['ASP_HW4_Problem_5_Real', '.fig']);

disp('(ii) Plotting imag part of x ...')
figure;
for ii = 1 : M
    subplot(M,1,ii)
    plot(1:L, imag(x_tilde_o(ii,:)))
    xlabel('Time index (n)')
    ylabel('Value')
    title(['Imag part of x_', int2str(ii), '(n | Y_n)'])
    xlim([1 L])
    grid on
end
saveas(gcf, ['figure/ASP_HW4_Problem_5_Imag', '.jpg']);
saveas(gcf, ['ASP_HW4_Problem_5_Imag', '.fig']);

disp('(iii) Plotting magnitude of x ...')
figure;
for ii = 1 : M
    subplot(M,1,ii)
    plot(1:L, abs(x_tilde_o(ii,:)))
    xlabel('Time index (n)')
    ylabel('Value')
    title(['Magnitude of x_', int2str(ii), '(n | Y_n)'])
    xlim([1 L])
    grid on
end
saveas(gcf, ['figure/ASP_HW4_Problem_5_Mag', '.jpg']);
saveas(gcf, ['ASP_HW4_Problem_5_Mag', '.fig']);

disp('(iv) Plotting unwrapped phase of x ...')
figure;
for ii = 1 : M
    subplot(M,1,ii)
    plot(1:L, unwrap(angle(x_tilde_o(ii,:))))
    xlabel('Time index (n)')
    ylabel('Radian')
    title(['Unwrapped phase of x_', int2str(ii), '(n | Y_n)'])
    xlim([1 L])
    grid on
end
saveas(gcf, ['figure/ASP_HW4_Problem_5_Phase', '.jpg']);
saveas(gcf, ['ASP_HW4_Problem_5_Phase', '.fig']);


disp('(c) Estimating the unwrapped phase of x ...')
w = zeros(M, 1);
u_tilde = 0.2;
eps = 0.001;
dis = 1;

for k = 1 : M
    phase_x = unwrap(angle(x_tilde_o(k,:)));
    res_phase = phase_x(2:L) - phase_x(1);
    
    for x = 1 : L-1
        y = w(k) * x;
        e = res_phase(x) - y;
        w(k) = w(k) + u_tilde / (x^2 + eps) * x * e;
    end
end

w_pi = w / pi;

for i = 1 : M
    fprintf(strcat("The slope of unwrapped phase of x", num2str(i), "(n) is: " ...
            ,num2str(round(w_pi(i), 3)), " pi \n"));
end



function [K_o, x_n1_tilde_o, x_tilde_o] = KalmanFilter(F_n1_n, C_n, Q1_n, Q2_n, K_i, x_tilde_i, y_tilde_i)
    function G = KalmanGain(F_n1_n, C_n, Q2_n, K_i)
        R = C_n * K_i * C_n' + Q2_n;
        G = F_n1_n * K_i * C_n' * inv(R);
    end
    function x_n1_tilde = OneStepPredictor(F_n1_n, C_n, x_tilde_i, y_tilde_i, G)
        alpha = y_tilde_i - C_n * x_tilde_i;
        x_n1_tilde = F_n1_n * x_tilde_i + G * alpha;
    end
    function K_o = RiccatiEquation(F_n1_n, C_n, Q1_n, G, K_i)
        Kn = K_i - inv(F_n1_n) * G * C_n * K_i;
        K_o = F_n1_n * Kn * F_n1_n' + Q1_n; 
    end
    G = KalmanGain(F_n1_n, C_n, Q2_n, K_i);
    K_o = RiccatiEquation(F_n1_n, C_n, Q1_n, G, K_i);
    x_n1_tilde_o = OneStepPredictor(F_n1_n, C_n, x_tilde_i, y_tilde_i, G);
    x_tilde_o = inv(F_n1_n) * x_n1_tilde_o;
end