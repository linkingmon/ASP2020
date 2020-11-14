function [a, P, kappa] = ASP_Levinson_Durbin(r)

% • a is a MATLAB cell array of size M. The entries in a contain the coefficients
% of the forward prediction error filter. More specifically, we have a{1} = a_1,
% a{2} = a_2, and a{M} = a_M.
% • P is an (M+1)-by-1 vector for the prediction errors. We have P = [P0; P1; P2; ... ; PM]^T .
% • kappa is an M-by-1 vector for the reflection coefficients. We have kappa =
% [kappa_1; kappa_2; ... ; kappa_M]^T

% Check if r is valid
if imag(r(1)) ~= 0
   error('ERROR: r(0) must be real!'); 
end

if ( r(1) < sqrt(max(conj(r).*r)) )
   error('ERROR: r(0) must be largest in magnitude!'); 
end

% Parameters
[rsize, ~] = size(r);
M = rsize - 1;
rB_M = flip(conj(r(2:end))); % RB_M = [r(-M) ... r(-1)]^T = [r(M) ... r(1)]^T*

a = {};
P = [];
kappa = [];

% Init
a{1} = 1;
P(1) = r(1);

% Levinson-Durbin Algorithm
for iter = 1:M
    delta(iter) = transpose(rB_M(end-iter+1:end)) * a{iter};
    kappa(iter) = -delta(iter) / P(iter);
    a{iter+1} = [a{iter}; 0] + kappa(iter) * [0; flip(conj(a{iter}))];
    P(iter+1) = P(iter) * (1 - kappa(iter)*conj(kappa(iter)) );
end

end