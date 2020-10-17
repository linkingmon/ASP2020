function J = ASP_Wiener_MSE(R, w, p, sd2)
    s_1 = size(R, 1);
    s_2 = size(R, 2);
    if ~(s_1 == s_2) || ~(s_1*s_2 == length(R(:)))
        error('The matrix R should be a square matrix');
    end
    if min(eig(R)) < 0
        error('The matrix R must be positive semi-definite');
    end
    if ~(sd2 > 0) || ~(imag(sd2) == 0)
        error('The variable sd2 must be a positive real number');
    end
    J = w'*R*w - w'*p - p'*w + sd2;
end