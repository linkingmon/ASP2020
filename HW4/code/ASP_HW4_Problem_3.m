mkdir figure
disp('Loading numerical values of matX ...')
load('ASP_HW4_Problem_3.mat')

[N, L] = size(matX);
disp(['Number of antenna is ', int2str(N)])
disp(['Number of measurements is ', int2str(L)])
R_hat = matX*matX' / L;
disp('(a) The sample covariance matrix is:')
disp(R_hat)

disp('(b) All eigenvalue of R_hat is:')
disp(eig(R_hat))

disp('(c) Calculating MVDR Spectrum ...')
P_MVDR = [];
theta_grid = [-90 : 0.01 : 90];
R_hat_inv = inv(R_hat);
for ii = 1 : numel(theta_grid)
    theta = theta_grid(ii);
    a = exp(1i*pi*sind(theta)*[0:N-1].');
    P_MVDR = [P_MVDR 1/real(a'*R_hat_inv*a)];
end
plot(theta_grid, P_MVDR)
xlabel('DOA (degree)')
ylabel('P_{MVDR}')
title('MVDR Spectrum')
saveas(gcf, ['figure/ASP_HW4_Problem_3_MVDR', '.jpg']);
saveas(gcf, ['ASP_HW4_Problem_3_MVDR', '.fig']);

[pks, locs] = findpeaks(P_MVDR, 'Npeaks', 6, 'SortStr', 'descend');
theta_pks = theta_grid(locs);
theta_pks =  theta_pks(logical(theta_pks<=10 & theta_pks>=0));
theta_s = theta_pks(1);

disp('(d) Estimate theta_s from MVDR Spectrum ...')
disp('Theta_s is:')
disp(theta_s)
a_s = exp(1i*pi*sind(theta_s)*[0:N-1].');

disp('(e) Plot Beam pattern of different weights ...')
B_uniform = [];
B_steering = [];
B_MVDR = [];
w_uniform = ones(N,1)/N;
w_steering = a_s/N;
w_MVDR = R_hat_inv*a_s/(a_s'*R_hat_inv*a_s);
for ii = 1 : numel(theta_grid)
    theta = theta_grid(ii);
    a = exp(1i*pi*sind(theta)*[0:N-1].');
    B_uniform =  [B_uniform abs(w_uniform'*a)];
    B_steering = [B_steering abs(w_steering'*a)];
    B_MVDR = [B_MVDR abs(w_MVDR'*a)];
end
plot(theta_grid, B_uniform)
hold on
plot(theta_grid, B_steering)
plot(theta_grid, B_MVDR)
xlabel('\theta (degree)')
ylabel('|B_{\theta}(\theta)|')
title('Beam pattern')
legend({'Uniform', 'Steering', 'MVDR'})
xlim([-90 90])
saveas(gcf, ['figure/ASP_HW4_Problem_3_Beampattern', '.jpg']);
saveas(gcf, ['ASP_HW4_Problem_3_Beampattern', '.fig']);

disp('(f) Plot Beamform Output of different weights ...')
y_uniform = [];
y_steering = [];
y_MVDR = [];
y_uniform = w_uniform'*matX;
y_steering = w_steering'*matX;
y_MVDR = w_MVDR'*matX;

figure;
subplot(2,1,1)
plot(1:L, real(y_uniform), 'LineWidth', 1)
hold on
plot(1:L, real(y_steering), 'LineWidth', 1)
plot(1:L, real(y_MVDR), 'LineWidth', 1)
xlabel('n')
ylabel('y(n)')
title('Beamform Output (real part)')
legend({'Uniform', 'Steering', 'MVDR'})

subplot(2,1,2)
plot(1:L, imag(y_uniform), 'LineWidth', 1)
hold on
plot(1:L, imag(y_steering), 'LineWidth', 1)
plot(1:L, imag(y_MVDR), 'LineWidth', 1)
xlabel('n')
ylabel('y(n)')
title('Beamform Output (imag part)')
legend({'Uniform', 'Steering', 'MVDR'})
saveas(gcf, ['figure/ASP_HW4_Problem_3_Output', '.jpg']);
saveas(gcf, ['ASP_HW4_Problem_3_Output', '.fig']);
