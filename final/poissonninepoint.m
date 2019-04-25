function w = poissonninepoint(a, b, h, f)
    x = a:h:b;
    y = a:h:b;    
    n = length(x) - 1;
    
    %construct coefficient matrix
    toe1 = toeplitz([-10/3, 2/3, zeros(1,n-3)]);
    toe2 = toeplitz([2/3, 1/6, zeros(1, n-3)]);
    super = diag(diag(eye(n-2)), 1);
    sub = diag(diag(eye(n-2)), -1);
    
    A = kron(eye(n-1), toe1) + kron(sub, toe2) + kron(super, toe2);
    
    [xx, yy] = meshgrid(x(2:end-1),y(2:end-1));

    flatxx = reshape(yy,[(n-1)*(n-1), 1]);
    flatyy = reshape(xx,[(n-1)*(n-1), 1]);
    
    [bx, by] = meshgrid(x, y);
    
    flatbx = reshape(by, [(n+1)*(n+1), 1]);
    flatby = reshape(bx, [(n+1)*(n+1), 1]);
    
    
    
    b = zeros((n-1)*(n-1), 1);
    for i =1:(n-1)*(n-1)
        b(i) = h^2*(1/12*f(flatbx(i+1), flatby(i)) + ...
            1/12*f(flatbx(i), flatby(i+1)) + ...
            2/3*f(flatbx(i+1), flatby(i+1)) +...
            1/12 * f(flatbx(i+2), flatby(i+2)) +...
            1/12 * f(flatbx(i+1), flatby(i+2)));
    end
    
    w = A\b;
    
    fullw = zeros(length(x), length(y));
    fullw([2:n], [2:n]) = reshape(w, n-1, n-1)';
    [xx, yy] = meshgrid(x, y);
    surf(xx, yy, fullw);
end