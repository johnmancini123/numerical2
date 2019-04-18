function output = midpoint(a, b, N, yo, f, realF)
h = (b-a)/N;
w = yo;
t = a;
output = zeros(3, N);
for i = 1: N
   w = w + h*(f(t+h/2, w + (h/2)*f(t, w)));
   t = t+h;
   output(1, i) = t;
   output(2, i) = w;
   output(3, i) = realF(t);
end
return;
end