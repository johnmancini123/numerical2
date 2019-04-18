function output = AdamsBashford2(a, b, N, yo, f, realF)
    output = zeros(3, N);
    %filling initial values
    mVals = yo;
    h = (b-a)/N;
    t = a + h;
    
    temp = rungeKutta4(a, b, N, yo, f, realF); 
    mVals = horzcat(mVals, temp(2, 1:end));%the vector that holds previous m values
    %filling in initial values for output
    output(1, 1) = t;
    output(2, 1) = mVals(2);
    output(3, 1) = realF(t);
    %done filling in init values for output
    for i = 2:N
       w = mVals(i) + (h/2)*(3*f(t, mVals(i)) - f(t-h, mVals(i-1)));
       mVals(i+1) = w;
       t = t+h;
       output(1, i) = t;
       output(2, i) = w;
       output(3, i) = realF(t);
    end
    return;
end
