%% ASP_HW1_Problem_5d
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
w0_real_n_sample = 101;
w0_imag_n_sample = 81;
w1 = 1 + 1i;
w2 = 0.5;
w0_real_sample = linspace(-5, 5, w0_real_n_sample);
w0_imag_sample = linspace(-4, 4, w0_imag_n_sample);


% Collect sample output data 
J_sample = [];
for ii = 1 : w0_real_n_sample
    for jj = 1 : w0_imag_n_sample
        w_sample = [w0_real_sample(ii)+1i*w0_imag_sample(jj) w1 w2].';
        J_sample(ii, jj) = ASP_Wiener_MSE(R, w_sample, p, sd2);
    end
end

% Plot figure
[w0_real_mesh, w0_imag_mesh] = meshgrid(w0_real_sample, w0_imag_sample);
figure;
surf(w0_real_mesh, w0_imag_mesh, real(J_sample.'));
colorbar;
title('MSE - w_0');
xlabel('Re(w_0)');
ylabel('Im(w_0)');
zlabel('MSE');

saveas(gcf, 'figure/ASP_HW1_Problem_5d.jpg');
