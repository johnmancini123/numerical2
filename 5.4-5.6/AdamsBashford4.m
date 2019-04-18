function output = AdamsBashford4(a, b, N, yo, f, realF)
    output = zeros(3, N);
    %filling initial values
    mVals = yo;
    h = (b-a)/N;    
    t = a;
    
    
    temp = rungeKutta4(a, b, N, yo, f, realF); 
    mVals = horzcat(mVals, temp(2, 1:end));%the vector that holds previous m values

    %filling in initial values for output
    for i = 1:3
       t = t+h;
       output(1, i) = t;
       output(2, i) = mVals(i+1);
       output(3, i) = realF(t);
    end
    
    %done filling in init values for output
    for i = 4:N
       w = mVals(i) + (h/24)*(55*f(t, mVals(i)) - 59*f(t-h, mVals(i-1)) +37*f(t-2*h, mVals(i-2)) - 9*f(t-3*h, mVals(i-3)));
       mVals(i+1) = w;
       t = t+h;
       output(1, i) = t;
       output(2, i) = w;
       output(3, i) = realF(t);
    end
    return;
end