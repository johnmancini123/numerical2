function W = wavefinitediff(L, tMax, h, k, alpha, f, t_left_bc, t_right_bc, g, n)
    m = L/h + 1;
    lambda = alpha*k/h;
    
    w = zeros(m, n);
  
    w(1,1) = 0;
    w(end, 1) = 0;
    for i =2: m-1
       w(i, 1) = f((i-1)*h);
       w(i, 2) = (1-lambda^2)*f((i-1)*h) + .5*lambda^2*(f((i-2)*h) + f(i*h)) + k*g((i-1)*h);
    end
    
    
    for j = 2 : n-1
        for i = 2 : m-1
            w(i, j+1) = 2*(1-lambda^2)*w(i, j) + lambda^2*(w(i+1, j)+w(i-1, j)) - w(i, j-1);
        end
    end
    
    W = w(1:end, end);
    
    x = 0:h:L;
    y = 0:k:tMax;
    [xx, yy] = meshgrid(y, x);
    surf(xx, yy, w);
    shg;
    xlabel('x')
    ylabel('t')
    return
    
    