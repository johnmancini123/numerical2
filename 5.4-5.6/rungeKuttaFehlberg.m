function output = rungeKuttaFehlberg(a, b, yo, hmin, hmax, tolerance, f, realF)
h = hmax;
w = yo;
t = a;
flag = 1;
i = 1;

while(flag)
    %calculating k values
    k1 = h*f(t, w);
    k2 = h*f(t + h/4, w + k1/4);
    k3 = h*f(t + (3*h)/8, w + (3/32)*k1 + (9/32)*k2);
    k4 = h*f(t + (12/13)*h, w + (1932/2197)*k1 - (7200/2197)*k2 + (7296/2197)*k3);
    k5 = h*f(t + h, w + (439/216)*k1 - 8*k2 + (3680/513)*k3 - (845/4104)*k4);
    k6 = h*f(t + h/2, w - (8/27)*k1 + 2*k2 - (3544/2565)*k3 + (1859/4104)*k4 - (11/40)*k5);
    %R = |ui+1 - wi+1|
    R = (1/h)*abs(k1/360 - (128/4275)*k3 - (2197/75240)*k4 + k5/50 + (2*k6)/55);
    %if absolute error is less than tolerance
    if (R <= tolerance)%approximation does not need modified step size
        t = t+h;
        output(4, i) = h;
        output(3, i) = realF(t);
        w = w + (25/216)*k1 + (1408/2565)*k3 + (2197/4104)*k4 - k5/5;
        output(1, i) = t;
        output(2, i) = w;
        
        i = i + 1;
    end
    %q = |(E*h/(2*R))^1/4|
    q = .84*(tolerance/R)^(.25);
    %if q is very low
    if q <= .1
        h = .1*h;
    %if q is very high
    elseif (q >= 4)
        h = 4*h;    
    %else re adjust step size
    else
        h = q*h;
    end
    %h is bigger than max h step size
    if h > hmax
        h = hmax;
    end
    %program is done
    if t >= b
        flag = 0;
    %if we've exceeded our maximum bound
    elseif (t+h) > b
        h = b-t;
    elseif h < hmin
            flag = 0;
            fprintf("minimum h exceeded\n");
    end
 
end
return;
end


