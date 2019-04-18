function [T Y SOLT] = ShootingMethodSecant(alpha, beta, a, b, t0, t1, tempYB, tol, str, soltF)
    FLAG = 1;
    MAXITER = 100;
    iter = 0;
    while FLAG
        [X, MAT] = ode45(str, [a b], [alpha;t1;0;1]); %solving an IVP with ode45
        %str is a variable that represents function name 'fname'. It is how
        %I pass function to ode45
        if abs(beta - MAT(end, 1)) > tol %recalculate t if we havent met tolerance
            temp = t1;
            t1 = t1 - ((MAT(end, 1) - beta)*(t1 - t0))/(MAT(end, 1)-tempYB); %secant method
            t0 = temp;
            tempYB = MAT(end, 1);
        else
           FLAG = 0; 
        end
        if iter == MAXITER
            fprintf("MAX ITERATIONS HIT\n");
            FLAG = 0;
        end
        iter = iter + 1;
    end
    fprintf("Number of iterations taken = %d\n", iter);
    T = X;
    Y = (MAT(1:end, 1));
    SOLT = soltF(T(1:end));
end