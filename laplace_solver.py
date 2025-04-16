import numpy as np
import matplotlib.pyplot as plt


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

    # Question 1
    exact_sol(50, 50)
