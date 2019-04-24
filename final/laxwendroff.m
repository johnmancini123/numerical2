function w = laxwendroff(alpha, beta, h, k, tMax, x_bc, A)
V = A*h/k;

n = (beta-alpha)/h;
m = tMax/k;
x = alpha:h:3;

%initializing starting vector according to problem
w0 = x_bc(x);
temp = zeros(1, ((beta-3)/h));
w0 = horzcat(w0, temp); %initial vector
X = alpha:h:beta;
Y = w0;
plot(X, Y);


w = zeros(1, n); %approximation vector

for j = 1: m
    for i = 1: n-1
        %lax wendroff
        w(i) = w0(i)*(.5*V*(1+V)) + (1-V^2)*w0(i+1) - .5*(V*(1-V))*w0(i+2);
    end
    %right endpoint
    w(i+1) = w0(i+1); %update right endpoint
    w0(1) = 0; %update left endpoint, the function is 0 in this case
    w0(2:end) = w; % update w0 to wj, so we can now calculate wj+1
end

X1 = alpha:h:beta;
Y1 = horzcat(0, w);
%plot(X, Y, 'r');
%plot(X1, Y1, 'g');
plot(X, Y, 'g', X1, Y1, 'r');
end