function w = cranknicolsonmethod(L, h, T, k, t_left_bc, t_right_bc, x_low_bc, diffusion)
X = 0:h:L;
Y = k:k:T;

fullu = zeros(length(Y)+1, length(X)); %matrix that will old all the points to plot 
a = t_left_bc(Y);
%had to brute force this because it would not allow me to do what I did
%with a
b = zeros(1, length(Y));
for i = 1:length(Y)
   b(i) = t_right_bc(Y(i)); 
end
fullu(1:end-1, 1) = a;
fullu(1:end-1, end) = b;
n = (L)/h; %number of iterations in space
m = T/k; %nubmer of iterations in time
x = h:h:L-h; % a matrix of the x values

lambda = (diffusion*(k/(h^2)));%calculating lambda

v1A = zeros(1, n-1) + (1+lambda); %making a vector for diagnolization
v2A = zeros(1, n-2) + (-lambda/2);
A = diag(v1A); %creating the (n-1)X(n-1) matrix A with v as diagno
temp = diag(v2A, 1);
A = A + temp;
temp = diag(v2A, -1);
A = A + temp;

v1B = zeros(1, n-1) + (1-lambda); %vector for diagnolization
B = diag(v1B); %B is now matrix with v2 as diagnol
v2B = zeros(1, n-2) + (lambda/2);
temp = diag(v2B, 1);
B = B + temp;
temp = diag(v2B, -1);
B = B + temp;
%adding upper and lower diagnol


t = k;

w_0 = x_low_bc(x); %initial vector of values
w_0 = transpose(w_0);
b = B*w_0;
b(1) = b(1) + t_left_bc(t);
b(end) = b(end) + t_right_bc(t);

fullu(end, 1:end) = x_low_bc(X);
w = zeros(n-1, 1); %the vector that holds the solutions

for j = 1: m
    w = A\b;
    fullu(end-j, 2:end-1) = w;
    t = t+k;
    b = B*w;
    b(1) = b(1) + t_left_bc(t);
    b(end) = b(end) + t_right_bc(t);
end
w = transpose(w);

Y = 0:k:T;
[xx, yy] = meshgrid(X, Y);
surf(xx, yy, fullu);
shg;
xlabel('x');
ylabel('y');

return;