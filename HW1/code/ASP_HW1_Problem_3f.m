%% ASP_HW1_Problem_3f
close all;
clear, clc

freq_sample = linspace(-2, 2, 1000000);

S_x = @(z) (0.199/0.99/(z-0.1) + 0.199*z/0.99/(1-0.1*z) + 3.97/0.99);
S_y = @(z) (7.939^2+7.5221^2+7.939*7.5221*z+7.939*7.5221/z) / 15.8801^2 * S_x(z);

S_y_sample = [];

for ii = 1 : numel(freq_sample)
    S_y_sample(ii) = S_y(exp(1i*2*pi*freq_sample(ii)));
end

figure;
plot(freq_sample, abs(S_y_sample), 'r-');
xlabel('Frequency (Hz)');
ylabel('Magnitude of PSD of y(n)');
title('PSD - frequency');
saveas(gcf, 'figure/ASP_HW1_Problem_3f-1.jpg');

figure;
plot(freq_sample, angle(S_y_sample), 'b-');
xlabel('Frequency (Hz)');
ylabel('Phase of PSD of y(n)');
title('PSD - frequency');
saveas(gcf, 'figure/ASP_HW1_Problem_3f-2.jpg');
