%%
%11.3 3c
fprintf("SECTION 11.3 NUMBER 3C\n");
p = @(x) (-(x+1));
q = @(x) (2);
r = @(x) ((1-x^2)*exp(-x));
a = 0;
b = 1;
h = .1;
y0 = -1;
y1 = 0;
[t, yApprox] = FiniteDifferenceLinear(p, q, r, a, b, h, y0, y1);
N = (b-a-h)/h;
for i = 1: N
   fprintf("i = %d t = %.2f yi = %.15f\n", i, t(i), yApprox(i)); 
end

fprintf("\n\n");

%%
%11.3 7
fprintf("SECTION 11.3 NUMBER 7\n");
p = @(x) (0);
q = @(x) (1000/(3*10^7*625));
r = @(x) (100*x*(x-120)/(2*3*10^7*625));
a = 0;
b = 120;
h = 6;
y0 = 0;
y1 = 0;
[t, yApprox] = FiniteDifferenceLinear(p, q, r, a, b, h, y0, y1);
actual_f = @(x) (7.7042537*10^4*exp(2.309401*10^(-4)*x) +  7.9207462*10^4*exp(-2.3094010*10^(-4)*x) -4.1666666*10^(-3)*(x - 120)*x - 1.5625*10^5);
N = (b-a-h)/h;
a_f = zeros(1, N);
fprintf("PART A\n");
for i = 1: N
    a_f(i) = actual_f(t(i));
   fprintf("i = %d t = %.2f w%d = %.15f y%d = %.15f abs value = %.15f\n", i, t(i), i, yApprox(i), i, a_f(i), abs(yApprox(i) - a_f(i))); 
   
end
max_error = max(abs(yApprox(1, 1:end) - a_f(1, 1:end)));
fprintf("\n\nPART B\n");
fprintf("max error = %.10f\n", max_error);
if (max_error > .02)
    fprintf("Max error not within .2 inches on interval\n");
else
    fprintf("Max error is within .2 inches on interval\n");
end

fprintf("\n\nPART C\n");
w_vec = zeros(1, 121);
for i = 0: 120
    w_vec(i+1) = actual_f(i);
end
max_w = max(w_vec);
fprintf("Max of real function = %.10f\n", max_w);
if (max_w < (1/300))
    fprintf("The actual function meets the standard\n");
else
    fprintf("The actual function does not meet the standard\n");
end

h = 1;
[t, yApprox] = FiniteDifferenceLinear(p, q, r, a, b, h, y0, y1);
max_approx = max(yApprox);
fprintf("Max of approximation function = %.10f\n", max_approx);

if (max_approx < (1/300))
    fprintf("The approximation function meets the standard\n");
else
    fprintf("The approximation function does not meet the standard\n");
end

fprintf("\n\n");


%%
%%12.1 3a
fprintf("SECTION 12.1 NUMBER 3A\n");
a = 0;  %left endpoint of x
b = 1;  %right endpoint of x
c = 0;  %left endpoint of y
d = 1;  %right endpoint of y
f = @(x,y) 0;  %right-hand side function

g1 = @(x,y) 0;        %u(x, c) BC
g2 = @(x,y) x; %u(x, d) BC
g3 = @(x,y) 0;        %u(a, y) BC
g4 = @(x,y) y;   %u(b, y) BC
%define step size
h = .2;



u = FiniteDifferenceElliptical(a, b, c, d, f, g1, g2, g3, g4, h, h);
realF = @(x, y) x*y;
x = a+h;
y = c+h;
r = (b-a-h)/h;
c = (d-c-h)/h;
for i = 1: r
    for j = 1: c
        fprintf("x = %.2f y = %.2f u = %.5f actual solt = %.5f\n", x, y,  u(round((i-1)*c + j)), realF(x, y));
        x = x + h;
    end
    x = a + h;
    y = y + h;
end

%%
%%12.1 3a
fprintf("SECTION 12.1 NUMBER 7\n");
a = 0;  %left endpoint of x
b = 6;  %right endpoint of x
c = 0;  %left endpoint of y
d = 5;  %right endpoint of y
f = @(x,y) -1.5/1.04;  %right-hand side function

g1 = @(x,y) x*(6-x);        %u(x, c) BC
g2 = @(x,y) 0; %u(x, d) BC
g3 = @(x,y) y*(5-y);        %u(a, y) BC
g4 = @(x,y) 0;   %u(b, y) BC
%define step size
h = .4;
k = 1/3;

u = FiniteDifferenceElliptical(a, b, c, d, f, g1, g2, g3, g4, h, k);
disp(length(u));
x = a+h;
y = c+k;
r = (b-a-h)/h;
c = (d-c-k)/k;
for i = 1: r
    for j = 1: c
        fprintf("%.15f x = %.2f y = %.2f\n", u(round((i-1)*c + j)), x, y);
        x = x + h;
    end
    x = a + h;
    y = y + k;
end