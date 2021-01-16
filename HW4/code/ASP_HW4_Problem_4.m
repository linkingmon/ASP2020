mkdir figure
disp('Loading numerical values of matX ...')
load('ASP_HW4_Problem_3.mat')

[N, L] = size(matX);
disp(['Number of antenna is ', int2str(N)])
disp(['Number of measurements is ', int2str(L)])
R_hat = matX*matX' / L;
R_hat_inv = inv(R_hat);
disp('(a) The sample covariance matrix is:')
disp(R_hat)

disp('(b) All eigenvalue of R_hat is:')
[v, d] = schur(R_hat); % R_hat = v*d*v'
disp(diag(d))

shift = 1;
Us = v(:,N-shift:N);
ds = d(N-shift:N,N-shift:N);

Us_1 = Us(1:N-1,:);
Us_2 = Us(2:N,:);

P_LS = pinv(Us_1)*Us_2;

s = eig(P_LS);
theta = asind(angle(s)/pi);
theta_1 = theta(1);
theta_2 = theta(2);

C = [exp(1i*pi*sind(theta_1)*[0:N-1].') exp(1i*pi*sind(theta_2)*[0:N-1].')];
g = [1; 0];
w_LCMV = R_hat_inv*C*inv(C'*R_hat_inv*C)*g;
y_LCMV = w_LCMV'*matX;
save('ASP_HW4_y_hat_t.mat', 'y_LCMV')

hold on
subplot(2,1,1)
plot(1:L, real(y_LCMV), 'LineWidth', 1)
xlabel('n')
ylabel('y(n)')
title('Beamform Output (real part)')
legend({'LCMV'})

subplot(2,1,2)
plot(1:L, imag(y_LCMV), 'LineWidth', 1)
xlabel('n')
ylabel('y(n)')
title('Beamform Output (imag part)')
legend({'LCMV'})
saveas(gcf, ['figure/ASP_HW4_Problem_4', '.jpg']);
