function w = backwarddifference(L, h, T, k, t_left_bc, t_right_bc, x_low_bc, diffusion);
n = (L)/h; %number of iterations in space
m = T/k; %nubmer of iterations in time
x = h:h:L-h; % a matrix of the x values

lambda = (diffusion*(k/(h^2)));%calculating lambda

v = zeros(1, n-1) + (1+2*lambda); %making a vector of 1-2*lambda
A = diag(v); %creating the (n-1)X(n-1) matrix A with 1-2*lambda in diagnol
%adding the lambda to above and below main diagnol
for i = 1:n-2
   A(i, i+1) = -lambda; 
end
for i = 2:n-1
    A(i, i-1) = -lambda;
end
w_0 = x_low_bc(x); %initial vector of values
w_0 = transpose(w_0);
w = zeros(1, n-1); %the vector that holds the solutions
w = transpose(w);
for j = 1: m
    w = A\w_0;
    w_0 = w;
end
w = transpose(w);
return;