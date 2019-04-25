function a = poissonfinitedifference(a, b, c, d, f, gxl, gxr, gyl, gyr, n, m, MAX, tol)
    h = (b-a)/n;%step 1
    k = (d-c)/m;
    
    x = (a+h):h:(b-h); %step 2
    y = (c+k):k:(d-k);%step 3
    
    w = zeros(n-1, m-1); %step 4
    
    lambda = (h/k)^2;%step 5
    mew = 2*(1+lambda);
    l = 1;
    
    NORM  = -1;
    while l <= MAX%step 6
        
        %step 7
        z = (-h^2*f(x(1), y(end)) + gyl(y(end)) + lambda*(gxr(d)) + lambda*w(1, end-1) + w(2, end))/mew;
        NORM = abs(z - w(1, end));
        w(1, end) = z;
        
        %step 8
        for i = 2: n-2
            z = (-h^2*f(x(i), y(end)) + lambda*(gxr(d)) + lambda*w(i-1, end) + w(i+1, end) + lambda*w(i, end-1))/mew;
            if abs(w(end, end) - z) > NORM:
                NORM = abs(w(end, end) - z);
            end
            w(end, end) = z;
        end
        
        %step 9
        z = (-h^2*f(x(end), y(end)) + gyr(y(end)) + lambda*gxr(x(end)) + w(end-1, end) + lambda*w(end, end-1))/mew;
        if abs(w(end, end) - z) > NORM:
            NORM = abs(w(end, end) - z);
        end
            w(end, end) = z;
            
        %step 10
        for j = m-2:-1:2
            
            %step 11
            z = (-h^2*f(x(j), y(j)) + gyl(y(j)) + lambda*w(1, j+1) + w(2, j) + lambda*w(i, j-1))/mew;
            if abs(w(1, j) - z) > NORM:
                NORM = abs(w(1, j) - z);
            end
            w(1, j) = z;
            
            %step 12
            for i =2:n-2
                z = (-h^2*f(x(i), y(j)) + w(i-1, j) + lambda*w(i, j+1) + w(i+1, j) + lambda*w(i, j-1))/mew;
                if abs(w(i,j) - z) > NORM:
                    NORM = abs(w(i, j) - z);
                end
                w(i, j) = z;
            end
            
            %step 13
            z = (-h^2*f(x(end), y(j)) + gyr(y(j)) + w(end-1, j) + lambda*w(end, j+1) + lambda*w(end, j))/mew;
            if abs(w(end, j) - z) > NORM:
                NORM = abs(w(end, j) - z);
            end
            w(end, j) = z;
        end
        
        %step 14
        
    end
    
end