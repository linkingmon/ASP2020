%% ASP_HW1_Problem_5c
clear, clc

% Setting parameters
R = [1.1 0.5 0.1 ;
     0.5 1.1 0.5 ;
     0.1 0.5 1.1 ] ;
p = [ 0.5 ;
     -0.4 ;
     -0.2] ;
sd2 = 1.0;

% Specify weight vector
n_sample = 101;
w0_real_sample = linspace(-5, 5, n_sample);
w0_imag = 1;
w1 = -0.5 + 1i;
w2 = -1;

% Collect sample output data 
J_sample = [];
for ii = 1 : n_sample
    w_sample = [w0_real_sample(ii)+1i*w0_imag w1 w2].';
    J_sample(ii) = ASP_Wiener_MSE(R, w_sample, p, sd2);
end

% Find minimum MSE
[min_J_val, min_J_idx] = min(real(J_sample));


% Plot
figure;
hold on;
plot(w0_real_sample, real(J_sample), 'b-.', 'LineWidth', 0.5);
plot(w0_real_sample(min_J_idx), min_J_val, 'r.', 'MarkerSize', 20);
text(-1.5, 1, strtrim(['Re(w_0) = ', num2str(w0_real_sample(min_J_idx)) ', min MSE = ', num2str(min_J_val)]), ...
    'VerticalAlignment','bottom','Color','red','FontSize',12);
title('MSE - Re(w_0)');
xlabel('Re(w_0)');
ylabel('MSE J');
grid on;
hold off;
legend('MSE','Minimum MSE', 'Location', 'northeast');

saveas(gcf, 'figure/ASP_HW1_Problem_5c.jpg');