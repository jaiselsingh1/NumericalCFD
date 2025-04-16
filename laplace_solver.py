import numpy as np
import matplotlib.pyplot as plt


def thomas_algorithm(p, q, r, s):

    n = len(s)
    rm = np.zeros(n)
    sm = np.zeros(n)
    x = np.zeros(n)

    rm[0] = r[0] / q[0]
    sm[0] = s[0] / q[0]

    # Direct sweep
    for i in range(1, n):
        qi = q[i] - p[i] * rm[i - 1]
        rm[i] = r[i] / qi
        sm[i] = (s[i] - p[i] * sm[i - 1]) / qi
    x[-1] = sm[-1]

    # Inverse sweep
    for i in range(n - 2, -1, -1):
        x[i] = sm[i] - rm[i] * x[i + 1]

    return x


def exact_sol(ix, jx):

    # Initializing arrays
    n = 50
    x = np.linspace(0.0, L, ix)
    y = np.linspace(0.0, W, jx)

    # Boundary Conditions
    phi = np.zeros((ix, jx))
    phi[:, -1] = np.ones(ix)
    
    # Exact solution
    for j in range(1, jx-1):
        for i in range(1, ix-1):
            for n in range(1, 51):
                phi[i, j] += 2/np.pi * ((-1)**(n+1) + 1)/n * np.sin(n * np.pi * x[i] / L) * np.sinh(n * np.pi * y[j] / L) / np.sinh(n * np.pi * W / L)
                
    # Contour Plot
    plt.figure()
    X, Y = np.meshgrid(x, y)
    plt.contourf(X, Y, phi.T)
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('y')

    # Plot at phi(x, W/2)
    plt.figure()
    idy = np.argmin(np.abs(y - W/2))
    plt.plot(x, phi[:, idy])
    plt.xlabel('x')
    plt.ylabel('phi(x, W/2)')

    plt.show()


if __name__ == "__main__":
    
    L = 1.0
    W = 1.0
    
    '''
    # Question 1
    exact_sol(50, 50)
    ''' 