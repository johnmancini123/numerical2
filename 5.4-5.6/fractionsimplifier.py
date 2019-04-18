def gcd(x, y):
    if y == 0:
        return x
    return gcd(y, x%y)


def simplify(n1, d1, n2, d2):
    n = (n1*d2) - (n2*d1)
    d = (d1*d2)
    div = gcd(d, n)
    n /= div
    d /= div
    print(n1/d1 - n2/d2)
    print(n/d)
    return(n, d)

if __name__ == '__main__':
    a = int(input('Enter numerator 1: '))
    b = int(input('Enter denomenator 1: '))
    c = int(input('Enter numerator 2: '))
    d = int(input('Enter denomenator 2: '))
    print(simplify(a, b, c, d))
    z = input()
