%%
%number 1
fprintf("Number 1\n");
n = 15;
m = 15;
a = 0;
b = 6;
c = 0;
d = 5;
fun = @(x, y) (-1.5/1.04);
gxl = @(x) (x*(6-x));
gxr = @(x) (0);
gyl = @(y) (5*(5-y));
gyr = @(y) (0);
MAX = 100;
tol = 10^-3;
w = poissonfinitedifference(a, b, c, d, fun, gxl, gxr, gyl, gyr, n, m, MAX, tol)

%%
%number 2 
fprintf("Number 2\n");
L = 1;
m = 100;
tMax = 3;
n = 300;
alpha = 1;
t_left = @(t) (0);
t_right = t_left;
g = @(x) (0);
w = wavefinitediff(L, m, tMax, n, alpha, @f, t_left, t_right, g);
h = .01;
x = 0:h:L;
for i = 1: length(w)
   fprintf("w(%.2f, %.2f) = %.10f\n", x(i), tMax, w(i)); 
end
%%
%number 3 approximations look right graph looks bad
L = 1;
m = 20;
n = 4800;
T = 6;
x_low = @(x) (1);

t_left = @(t) exp(-5*t);
t_right = @(t) (abs(cos(2*t))*exp(-(1/2)*t));
alpha = 1;

w = cranknicolson(L, m, T, n, alpha, t_left, t_right, x_low); %calculating approximation
h = .05;
x = 0:h:L;
%displaying results

fprintf("Crank nicolson method for problem 3 on final exam\n");
for i=1:length(x)
    fprintf("w(%.2f, %.2f) = %.10f\n", x(i), T, w(i));
end

fprintf("\n\n");



%%
%Number 5 
alpha = 0; 
beta = 10;
tMax = 5;
x_bc = @(x) (sin(2*pi*x));
h = .005;
k = .005;
A = 1;
w = laxwendroff(alpha, beta, h, k, tMax, x_bc, A);


%%

