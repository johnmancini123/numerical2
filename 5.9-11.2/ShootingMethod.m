function [T Y SOLT] = ShootingMethod(alpha, beta, a, b, tStart, tol, str, soltF)
    t = tStart;
    FLAG = 1;
    MAXITER = 100;
    iter = 0;
    while FLAG
        [X, MAT] = ode45(str, [a b], [alpha;t;0;1]); %solving an IVP with ode45
        if abs(beta - MAT(end, 1)) > tol %recalculate t
            t = t - (MAT(end, 1) - beta)/MAT(end, 3); %newton's method
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