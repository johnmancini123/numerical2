function zPrime = f3d(x, z)
    zPrime = zeros(1, 4);
    zPrime(1) = z(2);
    zPrime(2) = 1/2*(1 - z(2)^2 - z(1)*sin(x));
    zPrime(3) = z(4);
    zPrime(4) = -1/2*z(1)*cos(x)*z(3) -z(2)*z(4);
    zPrime = transpose(zPrime);
end
