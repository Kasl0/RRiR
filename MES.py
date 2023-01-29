from math import pi
import numpy as np
import matplotlib.pyplot as plt


# Oblicza funkcję testującą e dla zadanego n i x
def e(x, n):

    H = 3
    h = H/N

    if n < N:
        if (n-1)*h < x < n*h:
            return x/h - n + 1
        if n*h < x < (n+1)*h:
            return n + 1 - x/h

    return 0


# Oblicza pochodną funkcji testującej e dla zadanego n i x
def ed(x, n):
    
    H = 3
    h = H/N

    if n < N:
        if (n-1)*h < x < n*h:
            return 1/h
        if n*h < x < (n+1)*h:
            return -1/h

    return 0


# Oblicza B( e(n1), e(n2) ) czyli wyprowadzoną na kartce formułę
# Całka obliczana metodą trapezów
def B(n1, n2):

    integral = 0

    H = 3
    points = 100
    h = H/points

    for i in range(points):
        integral += h * (ed(h*i, n1) * ed(h*i, n2) + ed(h*(i+1), n1) * ed(h*(i+1), n2)) / 2

    return integral * (-1)


# Oblicza B( u, e(n) ) czyli składową potrzebną do wyliczenia L2(n)
# Całka obliczana metodą trapezów
def Bu(n):

    integral = 0

    H = 3
    points = 100
    h = H/points

    for i in range(points):
        integral += h * ((-1/3) * ed(h*i, n) + (-1/3) * ed(h*(i+1), n)) / 2

    return integral * (-1)


# Oblicza L(n) czyli wyprowadzoną na kartce formułę, która jest składową potrzebną do wyliczenia L2(n)
# Całka obliczana metodą trapezów
def L(n):

    integral = 0

    points = 100
    h = 1/points

    for i in range(points):
        integral += h * (e(1 + h*i, n) + e(1 + h*(i+1), n)) / 2

    G = 6.67259 / 10**11
    return 4 * pi * G * integral


# Oblicza L2 za pomocą zdefiniowanych wcześniej funkcji
def L2(n):
    return L(n) - Bu(n)


# Główna funkcja
# Generuje macierze, rozwiązuje układ równań A*X=C, wyznacza szukaną funkcję i rysuje wykres
def run(no_elements):

    # Ustawia globalną zmienną N jako parametr MES
    global N
    N = no_elements


    # Układ równań do rozwiązania: A * X = C

    #           [ B(1, 1)    B(2, 1)   ...   B(n−1, 1)]
    # Gdzie A = [ B(1, 2)    B(2, 2)   ...   B(n−1, 2)]
    #           [   ...        ...     ...      ...   ]
    #           [ B(1, n-1) B(2, n-1)  ... B(n−1, n-1)]

    # Gdzie C = [ L2(1) L2(2) ... L2(n-1)]

    # A X to szukana macierz postaci = [ w(1)  w(2) ... w(n-1)]

    n = N-1

    # Generuje macierz A
    A = np.empty((n,n))
    for i in range(n):
        for j in range(n):
            A[i,j] = B(j+1, i+1)

    # Generuje macierz C
    C = np.empty(n)
    for i in range(n):
        C[i] = L2(i)

    # Rozwiązuje układ równań, czyli oblicza macierz X
    X = np.linalg.solve(A, C)

    # Wyznacza szukaną funkcję na podstawie wartości X i wyznaczonej na kartce funkcji u
    def φ(x):

        u = 5 - x/3

        result = u
        
        for i in range(n):
            result += X[i] * e(x, i+1)

        return result

    # Rysuje wykres φ(x) w przedziale x ∈ [0, 3]
    no_points = 1000
    x = np.linspace(0, 3, no_points)

    y = np.empty(no_points)
    for i in range(no_points):
        y[i] = φ(x[i])

    plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("φ(x)")
    plt.title("Potencjał grawitacyjny φ")
    plt.show()


# Odpalenie głównej funkcji z arugmentem t (ilość elementów) jako parametrem uruchomieniowym aplikacji
t = 23
run(t)
