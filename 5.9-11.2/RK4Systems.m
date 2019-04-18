function [T, Z, Y] = RK4Systems(a, b, m, N, alpha, f, soltF)
    h = (b-a)/N; %step size
    t = a; %starting index
    w = zeros(1, m); %approximations
    
    T = zeros(1, N+1);%return T vector
    Z = zeros(2, N+1);%return vector for approximations
    Y = zeros(1, N+1);%return vec for real function values
  
    T(1) = t; 
    for j = 1:m
        w(j) = alpha(j);%adding initial value to both vectors
        Z(j, 1) = w(j);
    end
    Y(1) = alpha(1);
     
    k = zeros(4, m);
    for i = 1:N
           k(1, 1:m) = h*f(t, w); %getting k1 values for every function in the system of ODES
           k(2, 1:m) = h*f(t + h/2, w + (1/2)*k(1, :));%k2i vals for fi in the system of ODES
           k(3, 1:m) = h*f(t + h/2, w + (1/2)*k(2, :));%k3i vals for fi
           k(4, 1:m) = h*f(t + h, w + k(3, :));%k4i vals for fi
        for j = 1:m
            w(j) = w(j) + (k(1, j) + 2*k(2, j) + 2*k(3, j) + k(4, j))/6;%approximation for each fi in system of ODES
            Z(j, i+1) = w(j); %updating return vec
        end

        t = a + i*h;%updating t
        Y(i+1) = soltF(t); %real function val
        T(i+1) = t; %updating return vec
    end 
end