function w = thetamodified(L, h, T, k, t_left_bc, t_right_bc, x_low_bc, diffusion, r, p, c, theta);
n = (L)/h; %number of iterations in space
m = T/k; %nubmer of iterations in time
x = h:h:L-h; % a matrix of the x values

lambda = ((k/(diffusion*h^2)));%calculating lambda

v1 = zeros(1, n-1) + (1+2*lambda*theta); %making a vector for diagnolization
A = diag(v1); %creating the (n-1)X(n-1) matrix A with v as diagno
v2 = zeros(1, n-1) + (1-2*(1-theta)*lambda); %vector for diagnolization
B = diag(v2); %B is now matrix with v2 as diagnol
%adding upper and lower diagnol
for i = 1:n-2
   A(i, i+1) = -lambda*theta; 
   B(i, i+1) = (1-theta)*lambda;
end
for i = 2:n-1
    A(i, i-1) = -lambda*theta;
    B(i, i-1) = (1-theta)*lambda;
end
w_0 = x_low_bc(x); %initial vector of values
w_0 = transpose(w_0);

constant_vector = ones(n-1, 1)*(k*r/(p*c));
b = B*w_0 + constant_vector;
w = zeros(n-1, 1); %the vector that holds the solutions

for j = 1: m
    w = A\b;
    b = B*w + constant_vector;
end
w = transpose(w);
return;