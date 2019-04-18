%%5b
%for approximating f'(8.1)
x = [8.1 8.3 8.5 8.7];
fx = [16.94410 17.56492 18.19056 18.82091];
h = .2;
xo = 8.1;
approx = EndPoint(h, fx, x, xo, 1);
fprintf("The approximation of f'(%.1f) using the left threepoint formula is %.6f\n", xo, approx);

%%
%for approximating f'(8.3)
xo = 8.3;
approx = MidPoint(h, fx, x, xo);
fprintf("The approximation of f'(%.1f) using the mid threepoint formula is %.6f\n", xo, approx);


%%
%for approximating f'(8.5)
xo = 8.5;
approx = MidPoint(h, fx, x, xo);
fprintf("The approximation of f'(%.1f) using the mid threepoint formula is %.6f\n", xo, approx);

%%
%for approximating f'(8.7)
xo = 8.7;
approx = RightEndPoint(h, fx, x, xo, 4);
fprintf("The approximation of f'(%.1f) using the right threepoint formula is %.6f\n", xo, approx);

%%
%for approximating f''(1.3)
xo = 1.3;
h = .1;
fx = [11.59006 14.04276 16.86187];
approx = SecondDerivativeMidPoint(h, 2, fx);
fprintf("The approximation of f''(%.1f) using the midpoint formula with h = .1 is %.6f\n", xo, approx);

%%
%for approximating f''(1.3)
xo = 1.3;
h = .01;
fx = [13.78176 14.04276 14.30741];
approx = SecondDerivativeMidPoint(h, 2, fx);
fprintf("The approximation of f''(%.1f) using the midpoint formula with h = .01 is %.6f\n", xo, approx);
%%
function output = MidPoint(h, f, xVals, xo)
    %h is step size, f is function values at certain points
    %xVals is the x values of given data and xo is starting x point
    n = 0;
    for i = 1: 4
       if xo == xVals(i)
           n = i;
           break;
       end
    end
    output = (1/(2*h))*(f(n+1) - f(n-1));
    return;
end

function output = EndPoint(h, f, xVals, xo, n)
    %h is step size, f is function values at certain points

    %xVals is the x values of given data and xo is starting x point
    %n is index 
    output = (1/(2*h))*((-3*f(n)+ 4*(f(n+1)) - f(n+2)));
    return;
end

function output = RightEndPoint(h, f, xVals, xo, n)
    %h is step size, f is function values at certain points

    %xVals is the x values of given data and xo is starting x point
    %n is index 
    output = (1/(2*h))*((3*f(n) - 4*(f(n-1)) + f(n-2)));
    return;
end

function output = SecondDerivativeMidPoint(h, n, f)
    %h is step size, xo is value to approximate, n is starting index, and
    %f is f(x) values from given data
    output = (1/(h^2))*(f(n-1) -(2*(f(n))) + f(n+1));
    return;
end