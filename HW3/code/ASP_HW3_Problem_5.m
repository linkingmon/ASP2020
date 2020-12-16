mkdir figure;

plot_type = ''; % '', Jmin, excess, misadjust

Monte_Carlo_run(10, plot_type);
saveas(gcf, ['figure/ASP_HW3_Problem_5_R_10', plot_type, '.fig']);
Monte_Carlo_run(1000, plot_type);
saveas(gcf, ['figure/ASP_HW3_Problem_5_R_1000', plot_type, '.fig']);
close all;

function Monte_Carlo_run(realization, plot_type)
    load('ASP_HW3_Problem_5');
    
    [r, N] = size(matV);
    FIR_len = 50;
    h = @(n) ( (18/25)*(1/2).^n + (7/25)*(-1/3).^n );
    h = h([0:FIR_len-1]);
    M = 5;
    J_1 = zeros(N+FIR_len-1, 1);
    J_2 = zeros(N+FIR_len-1, 1);
    J_3 = zeros(N+FIR_len-1, 1);
    J_4 = zeros(N+FIR_len-1, 1);
    J_5 = zeros(N+FIR_len-1, 1);
    J_6 = zeros(N+FIR_len-1, 1);
    J_7 = zeros(N+FIR_len-1, 1);

    max_x_square = 0;
    for ii = 1 : realization
        v = matV(ii, :).';
        x = conv(h, v);
        max_x_square = max(max_x_square, max(x(1:500).*conj(x(1:500))));
        d = [v ; zeros(size(x,1)-size(v,1), 1)] ;
        J_1 = J_1 + ASP_LMS(x, d, M, 0.1)/realization;
        J_2 = J_2 + ASP_LMS(x, d, M, 0.2)/realization;
        J_3 = J_3 + ASP_NLMS(x, d, M, 0.2)/realization;
        J_4 = J_4 + ASP_NLMS(x, d, M, 0.8)/realization;
        J_5 = J_5 + ASP_RLS(x, d, M, 0.75, 0.01)/realization;
        J_6 = J_6 + ASP_RLS(x, d, M, 0.75, 0.1)/realization;
        J_7 = J_7 + ASP_RLS(x, d, M, 0.95, 0.01)/realization;
    end
    disp('Upper bound of u is')
    disp(1/max_x_square)
    J_min = 1.8835e-08;
    switch plot_type
    case 'excess'
        J_1(1:500) = J_1(1:500) - J_min;
        J_2(1:500) = J_2(1:500) - J_min;
        J_3(1:500) = J_3(1:500) - J_min;
        J_4(1:500) = J_4(1:500) - J_min;
        J_5(1:500) = J_5(1:500) - J_min;
        J_6(1:500) = J_6(1:500) - J_min;
        J_7(1:500) = J_7(1:500) - J_min;
    case 'misadjust'
        J_1(1:500) = (J_1(1:500) - J_min) / J_min;
        J_2(1:500) = (J_2(1:500) - J_min) / J_min;
        J_3(1:500) = (J_3(1:500) - J_min) / J_min;
        J_4(1:500) = (J_4(1:500) - J_min) / J_min;
        J_5(1:500) = (J_5(1:500) - J_min) / J_min;
        J_6(1:500) = (J_6(1:500) - J_min) / J_min;
        J_7(1:500) = (J_7(1:500) - J_min) / J_min;
    end
    figure
    semilogy([1:500], J_1(1:500), '-', 'Linewidth', 1.0);
    hold on
    semilogy([1:500], J_2(1:500), '-', 'Linewidth', 1.0);
    semilogy([1:500], J_3(1:500), '-', 'Linewidth', 1.0);
    semilogy([1:500], J_4(1:500), '-', 'Linewidth', 1.0);
    semilogy([1:500], J_5(1:500), '-', 'Linewidth', 1.0);
    semilogy([1:500], J_6(1:500), '-', 'Linewidth', 1.0);
    semilogy([1:500], J_7(1:500), '-', 'Linewidth', 1.0);
    if(strcmp(plot_type, 'Jmin'))
        semilogy([1:500], J_min*ones(500,1), '-', 'Linewidth', 1.0);
    end
    title('MSE');
    switch plot_type
    case 'excess'
        title('excess MSE');
    case 'misadjust'
        title('Misadjustment');
    end
    xlabel('n');
    ylabel('|e(n)|^2');
    legend('LMS u=0.1', 'LMS u=0.2', 'NLMS u=0.2', 'NLMS u=0.8', ...
        'RLS \lambda=0.75 \delta=0.01', 'RLS \lambda=0.75 \delta=0.1', 'RLS \lambda=0.95 \delta=0.01');
    if(strcmp(plot_type, 'Jmin'))
        legend('LMS u=0.1', 'LMS u=0.2', 'NLMS u=0.2', 'NLMS u=0.8', ...
        'RLS \lambda=0.75 \delta=0.01', 'RLS \lambda=0.75 \delta=0.1', 'RLS \lambda=0.95 \delta=0.01', 'J_{min}');
    end
    grid on
    hold off
end