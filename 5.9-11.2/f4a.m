function zPrime = f4a(x, z)
    zPrime = zeros(1, 4);
    zPrime(1) = z(2);
    zPrime(2) = z(1)^3 - z(1)*z(2);
    zPrime(3) = z(4);
    zPrime(4) = (3*z(1)^2 - z(2))*z(3) - z(1)*z(4);
    zPrime = transpose(zPrime);
end
