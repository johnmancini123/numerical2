function W = wavefinitediff(L, m, T, n, alpha, x_low, t_left, t_right, g)
    h = L/m;%step 1
    k = T/n;
    lambda = alpha*k/h;
    
    w = zeros(m+1, n+1); 
    w(1,1) = x_low(0);%step 3
    w(end, 1) = x_low(L);
    
    
    for i =2: m %step 4
       w(i, 1) = x_low((i-1)*h);
       w(i, 2) = (1-lambda^2)*x_low((i-1)*h) + .5*lambda^2*(x_low((i-2)*h) + x_low(i*h)) + k*g((i-1)*h);
    end
    
    
    for j = 2 : n %5
        for i = 2 : m
            w(i, j+1) = 2*(1-lambda^2)*w(i, j) + lambda^2*(w(i+1, j)+w(i-1, j)) - w(i, j-1);
        end
    end
    
    W = w(1:end, end);
    
    x = 0:h:L;
    y = 0:k:T;
    [xx, yy] = meshgrid(x, y);
    w = w';
    surf(xx, yy, w);
    shg;
    xlabel('x')
    ylabel('t')
    return
    
    