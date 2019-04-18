function output = rungeKutta4(a, b, N, yo, f, realF)
    h = (b-a)/N;
    w = yo;
    t = a;
    output = zeros(3, N);
    for i = 1:N
       k1 = f(t, w);
       halfH = h/2;
       k2 = f(t + halfH, w + halfH*k1);
       k3 = f(t + halfH, w + halfH*k2);
       k4 = f(t + h, w + h*k3);
       w = w + (halfH/3)*(k1 + 2*k2 + 2*k3 + k4);
       t = t+h;
       output(1, i) = t;
       output(2, i ) = w;
       output(3, i) = realF(t);
       
    end
    return;
end