from math import sqrt, pi
import numpy as np
import matplotlib.pyplot as plt

N = 23

def numericalIntegration(f):

    # points and weights read from Wikipedia:
    # https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Legendre_quadrature

    x_1 = 1 / sqrt(3)
    x_2 = -1 / sqrt(3)
    w_1 = 1
    w_2 = 1

    integral = w_1 * f(x_1) + w_2 * f(x_2)

    return integral


def e(x, n):

    H = 3
    h = H/N

    if n < N:
        if (n-1)*h < x < n*h:
            return x/h - n + 1
        if n*h < x < (n+1)*h:
            return n + 1 - x/h
    
    if n == N:
        if 0 < x < h:
            return 1 - x/h

    return 0


def ed(x, n):
    
    H = 3
    h = H/N

    if n < N:
        if (n-1)*h < x < n*h:
            return 1/h
        if n*h < x < (n+1)*h:
            return -1/h
    
    if n == N:
        if 0 < x < h:
            return -1/h

    return 0


def B(n1, n2):

    integral = 0

    H = 3
    points = 100
    h = H/points

    for i in range(points):
        integral += h * (ed(h*i, n1) * ed(h*i, n2) + ed(h*(i+1), n1) * ed(h*(i+1), n2)) / 2

    return integral * (-1)


def Bu(n):

    integral = 0

    H = 3
    points = 100
    h = H/points

    for i in range(points):
        integral += h * ((-1/3) * ed(h*i, n) + (-1/3) * ed(h*(i+1), n)) / 2

    return integral * (-1)


def L(n):

    integral = 0

    points = 100
    h = 1/points

    for i in range(points):
        integral += h * (e(1 + h*i, n) + e(1 + h*(i+1), n)) / 2

    G = 6.67259 / 10**11
    return 4 * pi * G * integral


def L2(n):
    return L(n) - Bu(n)


def run():

    n = N-1

    A = np.empty((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = B(j+1, i+1)

    C = np.empty(n)
    for i in range(n):
        C[i] = L2(i)

    X = np.linalg.solve(A, C)

    def u(x):

        result = 5 - x/3
        
        for i in range(n):
            result += X[i] * e(x, i+1)

        return result

    points = 1000

    x = np.linspace(0, 3, points)

    y = np.empty(points)
    for i in range(points):
        y[i] = u(x[i])

    plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("u(x)")
    plt.title("Potencjał grawitacyjny")
    plt.show()


run()
