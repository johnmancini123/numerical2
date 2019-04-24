%%
%Number 3 on final exam
diffusion = 1;
L = 1;
h = .05;
k = .00125;
tMax = 6;
x_low = @(x) (x-x+1);


t_right = @(t) (abs(cos(2*t))*exp(-1/2*t));
w = cranknicolsonmethod(L,h,tMax,k,t_left,t_right,x_low, 1); %calculating approximation
w = horzcat(t_left(T), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));

x = (0:h:L);
%displaying results

fprintf("Crank nicolson method for problem 3 on final exam\n");
for i=1:length(x)
    fprintf("w(%.2f, %.2f) = %.10f\n", x(i), T, w(i));
end

fprintf("\n\n");




%%
%Number 5 DONE
alpha = 0; 
beta = 10;
tMax = 5;
x_bc = @(x) (sin(2*pi*x));
h = .005;
k = .005;
A = 1;
w = laxwendroff(alpha, beta, h, k, tMax, x_bc, A);
%%
%number 2 CODE APPROXIMATES CORRECTLY, NEED TO GRAPH
L = 1;
tMax = 3;
h = .01;
k = .01;
alpha = 1;
t_left = @(t) (t-t+0);
t_right = t_left;
g = @(x) (x-x+0);
n = (tMax/k +1);
w = wavefinitediff(L, tMax, h, k, alpha,@f, t_left, t_right, g, n);

x = 0:h:L;
for i = 1: length(w)
   fprintf("w(%.2f, %.2f) = %.10f\n", x(i), tMax, w(i)); 
end

%%
L = 1;
tMax = .3;
h = .1;
k = .1;
alpha = 1;
t_left = @(t) (t-t+0);
t_right = t_left;
x_low = @(x) (sin(2*pi*x));
g = @(x) (2*pi*sin(2*pi*x));
n = int8(tMax/k +1);
w = wavefinitediff(L, tMax, h, k, alpha, x_low, t_left, t_right, g, n);

x = 0:h:L;
realf = @(x, t) (sin(2*pi*x)*(cos(2*pi*t) +sin(2*pi*t)));
for i = 1:length(x)
   fprintf("u(%.2f, %.2f) = %.10f w(%.2f, %.2f) = %.10f\n", x(i), tMax, realf(x(i), tMax), x(i), tMax, w(i));
end