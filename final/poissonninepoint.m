function w = poissonninepoint(a, b, h, n, f)
    %construct coefficient matrix
    A = zeros(n, n);
    
    v1 = zeros(1 ,n) + (-10/3);
    A = diag(v1);
    
    for i = 1: n-1
        if mod(i, 3) ~= 0
            A(i, i+1) = 2/3; 
            A(i+1, i) = 2/3;
        end
    end
    
    w = 1;
    v1 = zeros(1, n-3) + 2/3;
    temp = diag(v1, 3);
    A = A + temp;
    temp = diag(v1, -3);
    A = A + temp;
    
    for i = 1: n-3
        if mod(i, 3) ~= 0
            A(i+1, i+3) = 1/6; 
            A(i+3, i+1) = 1/6;
        end
    end
    
    
    disp(A);
end