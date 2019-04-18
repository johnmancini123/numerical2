function output = PredictorCorrector(a, b, N, yo, f, realF)
%4th order predictor corrector Adams method
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
        %this is how to be computationally efficient doing PC
        %here, ft0 means f(t-0*h, y_t-0*h)
       ft0 = f(t, mVals(i));
       ft1 = f(t-h, mVals(i-1));
       ft2 = f(t-2*h, mVals(i-2));
       wPredictor = mVals(i) + (h/24)*(55*ft0 - 59*ft1 +37*ft2 - 9*f(t-3*h, mVals(i-3)));
       t = t+h;
       w = mVals(i) + (h/24)*(9*f(t, wPredictor) + 19*ft0 - 5*ft1 + ft2);
       mVals(i+1) = w;
       output(1, i) = t;
       output(2, i) = w;
       output(3, i) = realF(t);
    end
    return;
end