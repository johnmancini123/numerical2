function u = FiniteDifferenceElliptical(a, b, c, d, f, g1, g2, g3, g4, h, k)

    %define vectors for grid
    x = a:h:b;
    y = c:k:d;

    %define number of subintervals
    n = length(x)-1;

    %create coefficient matrix for the (n-1)x(n-1) system 
    v = [-2 1 zeros(1, n-3)];
    D = sparse(toeplitz(v));
    A = kron(D, eye(n-1))+kron(eye(n-1), D);

    %define interior mesh points where u is unknown
    [xx, yy] = meshgrid(x(2:end-1),y(2:end-1));

    flatxx = reshape(yy,[(n-1)*(n-1), 1]);
    flatyy = reshape(xx,[(n-1)*(n-1), 1]);

    %create RHS vector for system
    b = zeros((n-1)*(n-1), 1);

    %update using RHS vector
    for i = 1:(n-1)*(n-1)
        b(i) = h^2*f(flatxx(i), flatyy(i));
    end

    %update using BCs at y(1) and y(end), g_1(x) and g_2(x)
    for i = 1:n-1
        b(i) = b(i) -g1(flatxx(i),y(1));
        b(end-i+1) = b(end-i+1)-g2(flatxx(end-i+1),y(end));
    end

    %update using BCs at x(1) and x(end), g_3(x) and g_4(x)
    for j = 1:n-1
        b(1+(j-1)*(n-1)) = b(1+(j-1)*(n-1))-g3(x(1),flatyy(1+(j-1)*(n-1)));
        b(n-1+(j-1)*(n-1)) = b(n-1+(j-1)*(n-1))-g4(x(end),flatyy(n-1+(j-1)*(n-1)));
    end

    %solve for u at interior mesh points
    u = A\b;

    %create mesh over the entire region for plotting
    [xx, yy] = meshgrid(x,y);

    %create a matrix to store u values on the entire region
    fullu = zeros(length(x), length(y));

    %update fullu with BC information
    for i = 1:n+1
        fullu(1, i) = g1(x(i),y(1));
        fullu(n+1, end-i+1) = g2(x(end-i+1),y(end));
        fullu(i, 1) = g3(x(1),y(i));
        fullu(i, n+1) = g4(x(end), y(i));
    end

    %update fullu with reshaped approximations at interior mesh points 
    fullu([2:n],[2:n]) = reshape(u, n-1, n-1)';

    %plot surface of approximations
    surf(xx,yy,fullu)
    shg
    xlabel('x')
    ylabel('y')
end