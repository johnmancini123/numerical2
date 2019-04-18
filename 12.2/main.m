%%
L = 2;
h = .4;
T = .5;
k = .1;
t_left = @(t) (0);
t_right = @(t) (0);
x_low = @(x) (sin(2*pi*x));
w = forwarddifference(L,h,T,k,t_left,t_right,x_low, 1); %calculating approximation
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
realF = @(x, t) (exp(-4*pi^2*t)*sin(2*pi*x));
x = 0:h:L;
u = realF(x, T); %creating vector for holding values of real functions at points xi and adding bc
n = int8((L+h)/h);

fprintf("Forward Difference method for Section 12.1 problem 5a with k = .1\n");
for i=1:n
    fprintf("u(%.2f, %.2f) = %.10f: w(%.2f, %.2f) = %.10f: abs_error = %.10f\n", x(i), T, u(i), x(i), T, w(i), abs(u(i)-w(i)));
end
fprintf("\n\n");

k = .05;
w = forwarddifference(L,h,T,k,t_left,t_right,x_low, 1); %calculating approximation
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
fprintf("Forward Difference method for Section 12.1 problem 5a with k = .05\n");
for i=1:n
    fprintf("u(%.2f, %.2f) = %.10f: w(%.2f, %.2f) = %.10f: abs_error = %.10f\n", x(i), T, u(i), x(i), T, w(i), abs(u(i)-w(i)));
end
fprintf("\n\n");
%%
L = 2;
h = .4;
T = .5;
k = .1;
t_left = @(t) (0);
t_right = @(t) (0);
x_low = @(x) (sin(2*pi*x));
w = backwarddifference(L,h,T,k,t_left,t_right,x_low, 1); %calculating approximation
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
realF = @(x, t) (exp(-4*pi^2*t)*sin(2*pi*x));
x = 0:h:L;
u = realF(x, T); %creating vector for holding values of real functions at points xi and adding bc
n = int8((L+h)/h);
fprintf("Backward Difference method for Section 12.1 problem 7a with k = .1\n");
for i=1:n
    fprintf("u(%.2f, %.2f) = %.10f: w(%.2f, %.2f) = %.10f: abs_error = %.10f\n", x(i), T, u(i), x(i), T, w(i), abs(u(i)-w(i)));
end
fprintf("\n\n");

k = .05;
w = backwarddifference(L,h,T,k,t_left,t_right,x_low, 1); %calculating approximation
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
fprintf("Backward Difference method for Section 12.1 problem 7a with k = .05\n");
for i=1:n
    fprintf("u(%.2f, %.2f) = %.10f: w(%.2f, %.2f) = %.10f: abs_error = %.10f\n", x(i), T, u(i), x(i), T, w(i), abs(u(i)-w(i)));
end

fprintf("\n\n");

%%
L = 2;
h = .4;
T = .5;
k = .1;
t_left = @(t) (0);
t_right = @(t) (0);
x_low = @(x) (sin(2*pi*x));
w = thetamethod(L,h,T,k,t_left,t_right,x_low, 1, .5); %calculating approximation
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
realF = @(x, t) (exp(-4*pi^2*t)*sin(2*pi*x));
x = 0:h:L;
u = realF(x, T); %creating vector for holding values of real functions at points xi and adding bc
n = int8((L+h)/h);
fprintf("Theta method(theta=.5) for Section 12.1 problem 9a with k = .1\n");
for i=1:n
    fprintf("u(%.2f, %.2f) = %.10f: w(%.2f, %.2f) = %.10f: abs_error = %.10f\n", x(i), T, u(i), x(i), T, w(i), abs(u(i)-w(i)));
end
fprintf("\n\n");

k = .05;
w = thetamethod(L,h,T,k,t_left,t_right,x_low, 1, .5); %calculating approximation
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
fprintf("Theta method(theta=.5) for Section 12.1 problem 9a with k = .05\n");
for i=1:n
    fprintf("u(%.2f, %.2f) = %.10f: w(%.2f, %.2f) = %.10f: abs_error = %.10f\n", x(i), T, u(i), x(i), T, w(i), abs(u(i)-w(i)));
end

fprintf("\n\n");

%%
L = 1.5;
diffusion = 1.04;
p = 10.6;
c = .056;
r = 5;
t_left = @(t) (0);
t_right = @(t) (0);
x_low = @(x) (sin(pi*x/1.5));
h = .15;
k = .0225;
T =  .54;
w = forwarddifferencemodified(L, h, T, k, t_left, t_right, x_low, diffusion, r, p, c);

n = (L+h)/h;
x = 0:h:L;
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
fprintf("Modifed forward difference with k = .0225 for question 13\n");
for i=1:n
    fprintf("w(%.2f, %.2f) = %.10f\n", x(i), T, w(i));
end

fprintf("\n\n");


%%
L = 1.5;
diffusion = 1.04;
p = 10.6;
c = .056;
r = 5;
t_left = @(t) (0);
t_right = @(t) (0);
x_low = @(x) (sin(pi*x/1.5));
h = .15;
k = .0225;
T = .54;
w = backwarddifferencemodified(L, h, T, k, t_left, t_right, x_low, diffusion, r, p, c);

n = (L+h)/h;
x = 0:h:L;
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
fprintf("Modifed backward difference with k = .0225 for question 13\n");
for i=1:n
    fprintf("w(%.2f, %.2f) = %.10f\n", x(i), T, w(i));
end

fprintf("\n\n");

%%
theta = 1/2;
L = 1.5;
diffusion = 1.04;
p = 10.6;
c = .056;
r = 5;
t_left = @(t) (0);
t_right = @(t) (0);
x_low = @(x) (sin(pi*x/1.5));
h = .15;
k = .0225;
T = .54;
w = thetamodified(L, h, T, k, t_left, t_right, x_low, diffusion, r, p, c, theta);

n = (L+h)/h;
x = 0:h:L;
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
fprintf("Modifed theta (theta=.5) with k = .0225 for question 13\n");
for i=1:n
    fprintf("w(%.2f, %.2f) = %.10f\n", x(i), T, w(i));
end

fprintf("\n\n");
