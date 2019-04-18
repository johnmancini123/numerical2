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
    fprintf("u(x%.2f, %.2f) = %.10f: w(x%.2f, %.2f) = %.10f: abs_error = %.10f\n", x(i), T, u(i), x(i), T, w(i), abs(u(i)-w(i)));
end
fprintf("\n\n");

k = .05;
w = forwarddifference(L,h,T,k,t_left,t_right,x_low, 1); %calculating approximation
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
fprintf("Forward Difference method for Section 12.1 problem 5a with k = .05\n");
for i=1:n
    fprintf("u(x%.2f, %.2f) = %.10f: w(x%.2f, %.2f) = %.10f: abs_error = %.10f\n", x(i), T, u(i), x(i), T, w(i), abs(u(i)-w(i)));
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
    fprintf("u(x%.2f, %.2f) = %.10f: w(x%.2f, %.2f) = %.10f: abs_error = %.10f\n", x(i), T, u(i), x(i), T, w(i), abs(u(i)-w(i)));
end
fprintf("\n\n");

k = .05;
w = backwarddifference(L,h,T,k,t_left,t_right,x_low, 1); %calculating approximation
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
fprintf("Backward Difference method for Section 12.1 problem 7a with k = .05\n");
for i=1:n
    fprintf("u(x%.2f, %.2f) = %.10f: w(x%.2f, %.2f) = %.10f: abs_error = %.10f\n", x(i), T, u(i), x(i), T, w(i), abs(u(i)-w(i)));
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
w = theta(L,h,T,k,t_left,t_right,x_low, 1, .5); %calculating approximation
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
realF = @(x, t) (exp(-4*pi^2*t)*sin(2*pi*x));
x = 0:h:L;
u = realF(x, T); %creating vector for holding values of real functions at points xi and adding bc
n = int8((L+h)/h);
fprintf("Theta method(theta=.5) for Section 12.1 problem 9a with k = .1\n");
for i=1:n
    fprintf("u(x%.2f, %.2f) = %.10f: w(x%.2f, %.2f) = %.10f: abs_error = %.10f\n", x(i), T, u(i), x(i), T, w(i), abs(u(i)-w(i)));
end
fprintf("\n\n");

k = .05;
w = theta(L,h,T,k,t_left,t_right,x_low, 1, .5); %calculating approximation
w = horzcat(t_left(0), w); %adding boundary conditions to the vector
w = horzcat(w, t_right(T));
fprintf("Theta method(theta=.5) for Section 12.1 problem 9a with k = .05\n");
for i=1:n
    fprintf("u(x%.2f, %.2f) = %.10f: w(x%.2f, %.2f) = %.10f: abs_error = %.10f\n", x(i), T, u(i), x(i), T, w(i), abs(u(i)-w(i)));
end

fprintf("\n\n");

%%
