%%
%Nothing
fprintf("John Mancini\nHomework 2\n 2/7/2019\n");
%%
%Section 5.2 
%5a and 7a
fty = @(t, y) (y/t - (y/t)^2);
a = 1;
b = 2;
N = 10;
init = 1;
yt = @(t) (t/(1+log(t)));


approx = Euler(fty, a,b, N, init, yt);

fprintf("Euler's method for 5a and 7a\n");

for i = 1: N
    fprintf("t = %.1f y = %.8f w = %.8f error = %.8f\n", approx(1, i), approx(2, i), approx(3, i), abs(approx(2,i)-approx(3,i)));
end

%%
%17
fty = @(t, y) (.1*.02*(1 - y));
a = 0;
b = 50;
N = 50;
init = .01;
yt = @(t) (1+(.01  - 1)*exp(1)^(-.1*.02*t));
approx = Euler(fty, a,b, N, init, yt);
fprintf("Euler's method for problem 17\n");
for i = 1: N
    fprintf("t = %.1f y = %.8f w = %.8f error = %.8f\n", approx(1, i), approx(2, i), approx(3, i), abs(approx(2,i)-approx(3,i)));
end


%%
%Section 5.3
%6b
fty = @(t,y) (y^2/(1+t));
a = 1;
b = 2;
N = 10;
init = -(log(2))^-1;
n = 2;
approx = Taylors(fty, a, b, N, init, n);
fprintf("Taylor method order 2 for problem 6b\n");
for i = 1: N
    fprintf("t = %.1f w = %.8f\n", approx(1, i), approx(2, i));
end
%8b
n = 4;
approx = Taylors(fty, a, b, N, init, n);
fprintf("\n\nTaylor method order 4 for problem 8b\n");
for i = 1: N
    fprintf("t = %.1f w = %.8f\n", approx(1, i), approx(2, i));
end
%%
%11
fty = @(t,y) (1 + t*sin(t*y));
a = 0;
b = 2;
N = 20;
init = 0;
n = 2;
approx = Taylors(fty, a, b, N, init, n);
fprintf("Taylors method of order 2 for problem 11\n");
for i = 1: N
    fprintf("t = %.1f w = %.8f\n", approx(1, i), approx(2, i));
end
