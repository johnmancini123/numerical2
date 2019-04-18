function output = rungeKuttaVerner(a, b, yo, hmin, hmax, tolerance, f, realF)
h = hmax;
w = yo;
t = a;
flag = 1;
i = 1;

while(flag)
    %calculating k values
    k1 = h*f(t, w);
    k2 = h*f(t + h/6, w + k1/6);
    k3 = h*f(t + (4*h)/15, w + (4/75)*k1 + (16/75)*k2);
    k4 = h*f(t + (2/3)*h, w + (5/6)*k1 - (8/3)*k2 + (5/2)*k3);
    k5 = h*f(t + (5/6)*h, w - (165/64)*k1 + (55/6)*k2 - (425/64)*k3 + (85/96)*k4);
    k6 = h*f(t + h, w + (12/5)*k1 - 8*k2 + (4015/612)*k3 - (11/36)*k4 + (88/255)*k5);
    k7 = h*f(t + h/15, w - (8263/15000)*k1 + (124/75)*k2 - (643/680)*k3 - (81/250)*k4 + (2484/10625)*k5);
    k8 = h*f(t + h, w + (3501/1720)*k1 - (300/43)*k2 + (297275/52632)*k3 - (319/2322)*k4 + (24068/84065)*k5 + (3850/26703)*k7);
    R = (1/h)*abs(-k1/160 - (125/17952)*k3 + k4/144 - (12/1955)*k5 - (3/44)*k6 + (125/11592)*k7 + (43/616)*k8);
        %if absolute error is less than tolerance
    if (R <= tolerance)%approximation does not need modified step size
        t = t+h;
        output(4, i) = h;
        output(3, i) = realF(t);
        w = w + (13/160)*k1 + (2375/5984)*k3 + (5/16)*k4 + (12/85)*k5 + (3/44)*k6;
        output(1, i) = t;
        output(2, i) = w;
        
        i = i + 1;
    end
    %q = |(E*h/(2*R))^1/4|
    q = .871*(tolerance/R)^(.2);
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
