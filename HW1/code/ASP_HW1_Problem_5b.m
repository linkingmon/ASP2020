%% ASP_HW1_Problem_5b
clear, clc

% Setting parameters
R = [1.1 0.5 0.1 ;
     0.5 1.1 0.5 ;
     0.1 0.5 1.1 ] ;
p = [ 0.5 ;
     -0.4 ;
     -0.2] ;
sd2 = 1.0;

% Calculate optimal weight vector
wopt = inv(R) * p;

% Calculate MSE J for optimal weight vector
J = ASP_Wiener_MSE(R, wopt, p, sd2);

str = 'The MSE of the optimal weight vector is: ';
Result = [str, num2str(J)];
disp(Result);
