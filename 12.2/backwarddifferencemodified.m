function w = backwarddifferencemodified(L, h, T, k, t_left_bc, t_right_bc, x_low_bc, diffusion, r, p, c);
n = (L)/h; %number of iterations in space
m = T/k; %nubmer of iterations in time
x = h:h:L-h; % a matrix of the x values

lambda = ((k/(diffusion*h^2)));%calculating lambda

v = zeros(1, n-1) + (1+2*lambda); %making a vector of 1-2*lambda
A = diag(v); %creating the (n-1)X(n-1) matrix A with 1-2*lambda in diagnol
%adding the lambda to above and below main diagnol
for i = 1:n-2
   A(i, i+1) = -lambda; 
end
for i = 2:n-1
    A(i, i-1) = -lambda;
end

w_0 = transpose(x_low_bc(x)); %initial vector of values

constant_vector = ones(n-1, 1)*(-k*r/(p*c)); %the vector that holds the constant
v = w_0-constant_vector;
w = zeros(n-1, 1); %the vector that holds the solutions
for j = 1: m
    w = A\v;
    v = w-constant_vector;
end
w = transpose(w);
return;