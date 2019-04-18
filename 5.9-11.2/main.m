
%% 
%5.9 4a
fprintf("SECTION 5.9 QUESTION 4a\n");
a = 0;
b = 1;
m = 2;
N = 10;
alpha = [2 2];
f = @(t, u) [(u(2)); (3*u(2) - 2*u(1) + 6*exp(-t))];
soltF = @(t) (2*exp(2*t) - exp(t) + exp(-t));
[t, z, y] = RK4Systems(a, b, m, N, alpha, f, soltF);
for i = 1: N+1
   fprintf("t = %.2f z(1, %d) = %.15f z(2, %d) = %.15f; solt(1, %d) = %.15f\n", t(i), i, z(1, i), i, z(2, i), i, y(i)); 
end
fprintf("\n\n");
%%
%11.1 3c
fprintf("SECTION 11.1 3C\n");
a = 0;
b = 1;
m = 2;
N = 10;
alpha1 = [-1 0];
f1 = @(t, u) [(u(2)); (-(t+1)*u(2) + 2*u(1) + (1-t^2)*exp(-t))];
alpha2 = [0 1];
f2 = @(t, u) [(u(2)); (-(t+1)*u(2) + 2*u(1))];
soltF = @(t) (1);
[t, z1, y] = RK4Systems(a, b, m, N, alpha1, f1, soltF); %using my RK4 method
[t, z2, y] = RK4Systems(a, b, m, N, alpha2, f2, soltF); 
BETA = 0;
final = z1(1, 1:end) + (BETA - z1(1, end))*z2(1, 1:end)/z2(1, end);
fprintf("Using my own implementation of rk4 for systems:\n");
for i = 1: N+1
   fprintf("t = %.2f y(1, %d) = %.15f\n", t(i), i, final(i)); 
end

fprintf("\n");
t = [0:.1:1];
[t1,u1]=ode45(@(t,u) f1(t,u),t,[-1 0],[]);
[t2,u2]=ode45(@(t,u) f2(t,u),t,[0 1],[]);

u1 = transpose(u1);
u2 = transpose(u2);
final = u1(1, 1:end) + (BETA - u1(1, end))*u2(1, 1:end)/u2(1, end);
fprintf("Using ode45\n");
for i = 1: N+1
   fprintf("t = %.2f y(1, %d) = %.15f\n", t(i), i, final(i)); 
end
fprintf("\n\n");
%%

%11.1 4b
fprintf("SECTION 11.1 4b\n");
a = 0;
b = pi/4;
m = 2;
N = 5;
alpha1 = [0 0];
f1 = @(t, u) [(u(2)); (-4*u(1) + cos(t))];
alpha2 = [0 1];
f2 = @(t, u) [(u(2)); (-4*u(1))];
soltF = @(t) (-1/3*cos(2*t) -sqrt(2)/6*sin(2*t)+1/3*cos(t));
[t, z1, y] = RK4Systems(a, b, m, N, alpha1, f1, soltF); %using my RK4 method
[t, z2, y] = RK4Systems(a, b, m, N, alpha2, f2, soltF); 
BETA = 0;
final = z1(1, 1:end) + (BETA - z1(1, end))*z2(1, 1:end)/z2(1, end);
fprintf("Using my own implementation of rk4 for systems:\n");
for i = 1: N+1
   fprintf("t = %.2f y(1, %d) = %.15f solt(%d) = %.15f\n", t(i), i, final(i), i, y(i)); 
end

fprintf("\n");
t = [0:pi/20:pi/4];
[t1,u1]=ode45(@(t,u) f1(t,u),t,[0 0],[]);
[t2,u2]=ode45(@(t,u) f2(t,u),t,[0 1],[]);

u1 = transpose(u1);
u2 = transpose(u2);
final = u1(1, 1:end) + (BETA - u1(1, end))*u2(1, 1:end)/u2(1, end);
fprintf("Using ode45\n");
for i = 1: N+1
   fprintf("t = %.2f y(1, %d) = %.15f\n", t(i), i, final(i)); 
end

fprintf("\n\n");

%%
%11.2 3d
fprintf("SECTION 11.2 3d\n");
alpha = 2;
beta = 2;
a = 0;
b = pi;
t = (beta-alpha)/(b-a);
soltF = @(t) (2 + sin(t));
tol = 10^-4;
fprintf("Using Shooting method with newton's method for root finding\n");
[T Y SOLT] = ShootingMethod(alpha, beta, a, b, t, tol, 'f3d', soltF);
for i = 1: length(T)
   fprintf("t = %.2f y(1, %d) = %.15f solt(1, %d) = %.15f\n", T(i), i, Y(i), i, SOLT(i)); 
end
fprintf("\n\n");
%%
%11.2 4a
fprintf("SECTION 11.2 4a\n");
alpha = .5;
beta = (1/3);
a = 1;
b = 2;
t = (beta-alpha)/(b-a);
soltF = @(t) ((t+1).^-1);
tol = 10^-4;
fprintf("Using Shooting method with newton's method for root finding\n");
[T Y SOLT] = ShootingMethod(alpha, beta, a, b, t, tol, 'f4a', soltF);
for i = 1: length(T)
   fprintf("t = %.2f y(1, %d) = %.15f solt(1, %d) = %.15f\n", T(i), i, Y(i), i, SOLT(i)); 
end
fprintf("\n\n");

%%
fprintf("SECTION 11.2 6b for 4a\n");
alpha = .5;
beta = (1/3);
a = 1;
b = 2;
t0 = (beta-alpha)/(b-a);
[X, MAT] = ode45('f4a', [a b], [alpha;t0;0;1]); %solving an IVP with ode45
tempYB = MAT(end, 1);
t1 = t0 + (beta - tempYB)/(b-a);

soltF = @(t) ((t+1).^-1);
tol = 10^-4;
fprintf("Using Shooting method with secant method for root finding\n");
[T Y SOLT] = ShootingMethodSecant(alpha, beta, a, b, t0, t1, tempYB, tol, 'f4a', soltF);
for i = 1: length(T)
   fprintf("t = %.2f y(1, %d) = %.15f solt(1, %d) = %.15f\n", T(i), i, Y(i), i, SOLT(i)); 
end
fprintf("\n\n");