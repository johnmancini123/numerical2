

%%
f = @(t, y) ((t+2*t^3)*y^3 - t*y);
a = 0;
b = 2;
hmax = .5;
hmin = .02;
tol = 10^(-6);
yo = 1/3;
realF = @(t) (3 + 2*t^2 + 6*exp(1)^(t^2))^(-.5);
[t_, w_ ,y_] = PredictorCorrectorAdaptive(a, b, yo, f, realF, tol, hmin, hmax);

j = 2;
while (t_(j) ~= 0)
   j = j +1; 
end

fprintf("Implementing adams method with variable step size for nubmer 3 on test\n");
for i = 1: j-1
   fprintf("i = %d t = %.10f w = %.10f error = %.10f\n", i, t_(i), w_(i), y_(i)); 
    
end