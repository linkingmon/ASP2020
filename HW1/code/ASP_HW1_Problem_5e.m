%% ASP_HW1_Problem_5e
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
w0_real_sample = linspace(-2, 2, n_sample);
w0_imag = 0;
w1 = -0.7683;
w2_imag = 0;
w2_real_sample = linspace(-2, 2, n_sample);

J_sample = [];
for ii = 1 : n_sample
    for jj = 1 : n_sample
        w_sample = [w0_real_sample(ii)+1i*w0_imag w1 w2_real_sample(jj)+1i*w2_imag].';
        J_sample(ii, jj) = ASP_Wiener_MSE(R, w_sample, p, sd2);
    end
end

% Find minimum
min_J_val = min(J_sample(:));
[row_idx, col_idx] = find(J_sample == min_J_val);

% Plot figure
[w0_real_mesh, w2_real_mesh] = meshgrid(w0_real_sample, w2_real_sample);
figure;
hold on;
lev = [0.35, 0.6, 1, 2, 3, 4, 5];
contour(w0_real_mesh, w2_real_mesh, J_sample.', lev, 'ShowText','on');
plot(w0_real_sample(row_idx),w2_real_sample(col_idx),'r.', 'MarkerSize', 20);
mintxt = sprintf(['w_0 = ', num2str(real(w0_real_sample(row_idx))), ', w_2 = ',  num2str(real(w2_real_sample(col_idx))), ...
    ', \nmin MSE = ', num2str(min_J_val)]);
strValues = strtrim(mintxt);
text(0.6, 0.3, strValues,'VerticalAlignment','bottom','Color','red','FontSize',10);
title('MSE - Re(w_0) x Re(w_2)');
xlabel('Re(w_0)');
ylabel('Re(w_2)');
legend('MSE','Minimum MSE', 'Location', 'southeast');
hold off;

saveas(gcf, 'figure/ASP_HW1_Problem_5e.jpg');
