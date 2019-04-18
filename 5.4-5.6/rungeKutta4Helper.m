function [T, Y] = rungeKutta4Helper(a, h, yo, f)
%an rk4 method that only does 3 calculations as a helper for the predictor
%corrector adaptive method
    T = zeros(1, 4);
    Y = zeros(1, 4);

    w = yo;
    t = a;
    T(1) = t;
    Y(1) = w;
    for i = 2:4
       k1 = f(t, w);
       halfH = h/2;
       k2 = f(t + halfH, w + halfH*k1);
       k3 = f(t + halfH, w + halfH*k2);
       k4 = f(t + h, w + h*k3);
       w = w + (halfH/3)*(k1 + 2*k2 + 2*k3 + k4);
       t = t+h;
       T(i) = t;
       Y(i) = w;
       
    end
    return;
end