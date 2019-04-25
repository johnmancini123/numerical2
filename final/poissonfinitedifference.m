function a = poissonfinitedifference(a, b, c, d, f, gxl, gxr, gyl, gyr, h, k)
    x = a:h:b; %length are the same
    y = c:k:d;
    lambda = (h/k)^2;
    %defining the number of sub intervals
    n = length(x)-1; %number of sub intervals+1 in both x and y directionss
    
    b = zeros((n-1)*(n-1), 1);
    
    %define interior mesh points where u is unknown
    [xx, yy] = meshgrid(x(2:end-1),y(2:end-1));
    
    flatxx = reshape(yy,[(n-1)*(n-1), 1]);
    flatyy = reshape(xx,[(n-1)*(n-1), 1]);

    %update using RHS vector
    for i = 1:(n-1)*(n-1)
        b(i) = h^2*f(flatxx(i), flatyy(i));
    end
    
    
    
    %TODO create a coefficient matrix of system
end