function w = poissonfinitedifference(a, b, c, d, f, gxl, gxr, gyl, gyr, n, m, MAX, TOL)
    h = (b-a)/n;%step 1
    k = (d-c)/m;
    x = (a+h):h:(b-h); %step 2
    y = (c+k):k:(d-k);%step 3
    
    w = zeros(n-1, m-1); %step 4
    
    lambda = (h/k)^2;%step 5
    mew = 2*(1+lambda);
    l = 1;
    
   
    while l <= MAX%step 6
        
        %step 7
        z = (-h^2*f(x(1), y(end)) + gyl(y(end)) + lambda*(gxr(x(1))) + lambda*w(1, end-1) + w(2, end))/mew;
        NORM = abs(z - w(1, end));
        w(1, end) = z;
        
        %step 8
        for i = 2: n-2
            z = (-h^2*f(x(i), y(end)) + lambda*(gxr(x(i))) + w(i-1, end) + w(i+1, end) + lambda*w(i, end-1))/mew;
            if abs(w(end, end) - z) > NORM
                NORM = abs(w(i, end) - z);
            end
            w(i, end) = z;
        end
        
        %step 9
        z = (-h^2*f(x(end), y(end)) + gyr(y(end)) + lambda*gxr(x(end)) + w(end-1, end) + lambda*w(end, end-1))/mew;
        if abs(w(end, end) - z) > NORM
            NORM = abs(w(end, end) - z);
        end
            w(end, end) = z;
            
        %step 10
        for j = m-2:-1:2
            
            %step 11
            z = (-h^2*f(x(1), y(j)) + gyl(y(j)) + lambda*w(1, j+1) + w(2, j) + lambda*w(1, j-1))/mew;
            if abs(w(1, j) - z) > NORM
                NORM = abs(w(1, j) - z);
            end
            w(1, j) = z;
            
            %step 12
            for i =2:n-2
                z = (-h^2*f(x(i), y(j)) + w(i-1, j) + lambda*w(i, j+1) + w(i+1, j) + lambda*w(i, j-1))/mew;
                if abs(w(i,j) - z) > NORM
                    NORM = abs(w(i, j) - z);
                end
                w(i, j) = z;
            end
            
            %step 13
            z = (-h^2*f(x(end), y(j)) + gyr(y(j)) + w(end-1, j) + lambda*w(end, j+1) + lambda*w(end, j-1))/mew;
            if abs(w(end, j) - z) > NORM
                NORM = abs(w(end, j) - z);
            end
            w(end, j) = z;
        end
        
        %step 14
        z = (-h^2*f(x(1), y(1)) + gyl(y(1)) + lambda*gxr(x(1)) + lambda*w(1, 2) + w(2, 1))/mew;
            if abs(w(1, 1) - z) > NORM
                NORM = abs(w(1, 1) - z);
            end
        w(1, 1) = z;
        
        %step 15
        for i = 2:n-2
            z = (-h^2*f(x(i), y(1)) + lambda*gxr(x(i)) + w(i-1, 1) + lambda*w(i, 2) + w(i+1, 1))/mew;
            if abs(w(i, 1) - z) > NORM
                NORM = abs(w(i, 1) - z);
            end
            w(i, 1) = z;
        end
        
        %step 16
        z = (-h^2*f(x(end), y(1)) + gyr(y(1)) + lambda*gxr(x(end)) + w(end-1, 1) + lambda*w(end, 2))/mew;
            if abs(w(end, 1) - z) > NORM
                NORM = abs(w(end, 1) - z);
            end
        w(end, 1) = z;
        
        %step 17
        if NORM <= TOL
           %step 18 output results 
           return; 
        end
        
        l = l+1;
    end
    
    if (l>=MAX)
        fprintf("MAX ITER REACHED\n");
    end
    
    
end