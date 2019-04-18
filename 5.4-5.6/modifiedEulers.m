function output = modifiedEulers(a, b, N, yo, f, realF)
%a is low b is high N is number of steps
%yo is initial value and f is function = y'(t) 
%and realF is the function we're trying to approximate
output = zeros(3, N);
h = (b-a)/N;
w = yo;
t = a;
for i = 1: N
   w = w + (h/2)*(f(t, w) + f(t+h, w + h*(f(t,w))));
   t = t+h;
   output(1, i) = t;
   output(2, i) = w;   
   output(3, i) = realF(t);
end
return;
end