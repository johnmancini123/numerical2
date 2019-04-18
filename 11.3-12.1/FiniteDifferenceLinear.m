function [t_retVec, Y] = FiniteDifferenceLinear(p, q, r, a, b, h, alpha, beta)
    N = (b-a-h)/h; %getting N
    T = zeros(N, N); %T matrice 
    B = zeros(1, N); %b value matrice
    t_retVec = zeros(1, N); % returning the vector of t values
    t_retVec(1) = a+h;
    t = a+h; %starting t is a + h
    T(1, 1) = -(2+h^2*q(t)); %first row of matrix
    T(1, 2) = (1 - h*p(t)/2);
    B(1) = (h^2*r(t) - (1+h*p(t)/2)*alpha); %first input to b vector
    t_retVec(N) = b-h;
    t = b - h; %final t value in between
    T(N, N-1) = (1 + h*p(t)/2); %last row of T matrix
    T(N, N) = -(2 + h^2*q(t));
    B(N) = h^2*r(t) - (1 - h*p(t)/2)*beta; % last input to b vector
  
    t = a + 2*h; %starting value of t
    for i = 2: N-1
            T(i, i-1) = 1 + h*p(t)/2;%calculating 3 values per row
            T(i, i) = -(2 + h^2*q(t)); 
            T(i, i+1) = 1 - h*p(t)/2;
            B(i) = h^2*r(t);%calculating b value in column
            t_retVec(i) = t; %updating t vector
            t = t+h;%updating t
    end
    Y = linsolve(T, transpose(B)); %solving the system
end