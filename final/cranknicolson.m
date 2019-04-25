function w = cranknicolson(L, m, T, n, alpha, left_bc, right_bc, low_bc)
    fullu = zeros(n+1, m+1);

    h = L/m;%step 1
    k = T/n;
    N = n;
    lambda = alpha^2*k/h^2;
    w = zeros(1, m+1);
    w(1) = left_bc(0);
    w(end) = right_bc(0);

    for i = 2:m %step 2
        w(i) = low_bc((i-1)*h);
    end

    fullu(end, 1:end) = w;

    l = zeros(1, m-1);
    u = zeros(1, m-1);

    l(1) = 1+lambda; %step 3
    u(1) = -lambda/(2*l(1));

    for i = 2:m-2 %step 4
       l(i) = 1 + lambda + lambda*u(i-1)/2;
       u(i) = -lambda/(2*l(i));
    end

    l(end) = 1+lambda + lambda*u(end-1)/2; %step 5

    z = zeros(1, m-1);
    
    for j = 1:N %step 6
       t = j*k;
       z(1) = ((1-lambda)*w(2) + lambda/2*w(3) + lambda/2*w(1) + lambda/2*left_bc(t))/l(1); %step 7
       w(1) = left_bc(t);

       for i = 2:m-2  %step 8
           z(i) = ((1-lambda)*w(i+1) +lambda/2*(w(i+2) + w(i) + z(i-1)))/l(i);
       end
       
       z(end) = ((1-lambda)*w(end-1) + lambda/2*(w(end) + w(end-2) + z(end-1) + right_bc(t)))/l(end);
       w(end) = right_bc(t);



       w(end-1) = z(end); %step 9

       for i = m-1:-1:2 %step 10
           w(i) = z(i-1) - u(i-1)*w(i+1);
       end
 
       fullu(end-j, 1:end) = w;
    end
    
    
    x = 0:h:L;
    y = 0:k:T;
    [xx, yy] = meshgrid(x, y);
    
    mesh(xx, yy, fullu);
    shg;
    xlabel('x');
    ylabel('t');
end