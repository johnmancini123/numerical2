%%
%Implementing Euler's method
function output = Euler(f, a, b, N, init, actualFunc)
%Implementation of Euler's method where 
%f = f(t, y), a = left endpoint, b = right endpoint,
%N = number of steps, and init = value of initial condition
h = (b-a)/N;%calculating h
w = init;%setting w to the initila value
t = a;%t should start at a

toRet = zeros(3, N);
for i = 1: N
   w = w + h*f(t, w); %approximate value from Euler's
   t = a + i*h; %updating t
   y = actualFunc(t); %the actual value at t
   toRet(1, i) = t;
   toRet(2, i) = y;
   toRet(3, i) = w;
end
output = toRet;
return;
end