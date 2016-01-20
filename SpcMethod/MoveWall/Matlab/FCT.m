%% Fast Chebyshev Transform
%% Introduction
%   This function is used to transform a series of data from
%   physical space into spectral space (or inverse), which is
%   useful in solving boundary problems by spectral method.
%% Theory
%   Suppose we have a series of data, {f_k}, k = 0,1,2,...,N,
%   which is generated from a function f(x) on Chebyshev-
%   Gauss-Lobatto Nodes x=cos(k*pi/N). Let T_k(x) be the k-th
%   Chebyshev polynomial, which is cos(k * arccos(x)). Now
%   we can expand f(x) as Chebyshev series,
%
%           f(x) = sum a_m * T_m(x), m = 0:N
%
%   Since T_m(x_k) = cos(m*n*pi/N), we can get
%
%           f_k = sum a_m * cos(m*n*pi/N), m = 0:N
%
%   Observing this formula, we can discover its relationship
%   with DFT. To make this clear, we construct a new sequence,
%   which is {f_0, f_1, ... , f_N , f_N-1, ... , f_2}, and perform
%   2N IDFT on it, we get
%
%           f_k = sum a_m * exp(2*pi*i * m*k/2*N), m = 0:2N
%
%   Since f_k = f_2N-k, we get
%
%           sum a_m * i * sin(2*pi * m*k/2*N) = 0, m = 0:2N
%
%   This indicates that the IDFT above only contains real part.
%   And we get
%
%           f_k = sum a_m * cos(pi*m*n/N), m = 0:2N
%
%   Note that exp(2*pi*i * m*(2N-k)/2N) = exp(2*pi*i * (2N-m)*k/2N),
%   we can get a_m = a_2N-m. So that
%
%           f_k = a_0 + sum 2*a_m*cos(m*k*pi/N) + a_N*cos(N*k*pi/N),
%                       m = 1:N-1
%   and we now get the answer for the formal problem.

%% Parameters
%   A ---------------- Input data in coloum
%   B ---------------- Output data in coloum
%   direction -------- 1 for physical -> spectral
%                     -1 for spectral -> physical

%% Author
%   Written by Luohao Wang on 2015-9-10
%   Contact : lh-wang13@mails.tsinghua.edu.cn

%% Acknowledgement
%   Inspired by fcgltran by Greg von Winckel, fixed a few bugs.
%   Website : http://www.mathworks.com/matlabcentral/fileexchange/4591-fast-chebyshev-transform--1d-

%% Code
function B = FCT(A, direction)

[N,~] = size(A);
if N == 1
    B = transpose(FCT(transpose(A), direction));
else
    if direction == 1 % physical -> spectral
        F=ifft([A(1:N,:);A(N-1:-1:2,:)]);
        B=([F(1,:); 2*F(2:(N-1),:); F(N,:)]);
    elseif direction == -1            % physical -> spectral
        F=fft([A(1,:);A(2:N-1,:)/2;A(N,:);A(N-1:-1:2,:)/2]);
        B=(F(1:N,:));
    else
        disp('Direction identifier can only be -1 or 1! Please check your program.');
    end
end
end