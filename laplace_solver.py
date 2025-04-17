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

    # Constants
    n = 50

    # Initializing arrays
    x = np.linspace(0.0, L, ix)
    y = np.linspace(0.0, W, jx)
    phi = np.zeros((ix, jx))

    # Boundary Conditions
    phi[:, -1] = np.ones(ix)
    
    # Exact solution
    for j in range(1, jx-1):
        for i in range(1, ix-1):
            for n in range(1, 51):
                phi[i, j] += 2/np.pi * ((-1)**(n+1) + 1)/n * np.sin(n * np.pi * x[i] / L) * np.sinh(n * np.pi * y[j] / L) / np.sinh(n * np.pi * W / L)
                
    # Contour Plot
    plt.figure('Exact Solution Contour Plot')
    X, Y = np.meshgrid(x, y)
    plt.contourf(X, Y, phi.T)
    plt.colorbar()
    plt.xlabel('x')
    plt.ylabel('y')

    # Plot at phi(x, W/2)
    plt.figure('Exact Solution Line Plot')
    idy = np.argmin(np.abs(y - W/2))
    plt.plot(x, phi[:, idy])
    plt.xlabel('x')
    plt.ylabel('phi(x, W/2)')

    return phi[:, idy]


def jacobi(ix, jx):

    # Constants
    res = 1e-5
    dx = L/(ix-1)
    dy = W/(jx-1)
    K = 2*(1/dx**2 + 1/dy**2) # 1/h^2
    max_iter = 10000

    # Initializing arrays
    x = np.linspace(0.0, L, ix)
    y = np.linspace(0.0, W, jx)
    phi = np.zeros((ix, jx))
    phi_new = np.zeros((ix, jx))
    max_res_list = []

    # Boundary Condition
    phi[:, -1] = np.ones(ix)
    phi_new[:, -1] = np.ones(ix)

    for k in range(max_iter):

        # Jacobi
        for j in range(1, jx-1):
            for i in range(1, ix-1):
                 phi_new[i, j] = 1/K * ((phi[i+1, j] + phi[i-1, j])/dx**2 + (phi[i, j+1] + phi[i, j-1])/dy**2)

        # Computing max residual
        max_res = np.max(np.abs(phi_new - phi))
        max_res_list.append(max_res)

        if max_res < res:
            break

        # Shifting computation window back
        phi = phi_new.copy()

    # Plotting
    plt.figure('Jacobi Max Residual')
    plt.semilogy(np.arange(k+1), max_res_list, label=f"ix = jx = {ix}")
    plt.xlabel('Iterations')
    plt.ylabel('Max Residual')
    plt.legend()
    
    plt.figure('Jacobi Line Plot')
    idy = np.argmin(np.abs(y - W/2))
    plt.plot(x, phi[:, idy], label=f"ix = jx = {ix}")
    plt.xlabel('x')
    plt.ylabel('phi(x, W/2)')
    if ix == 50:
        plt.plot(x, phi_exact, label="Exact Solution")
    plt.legend()


def gauss_seidel(ix, jx):

    # Constants
    res = 1e-5
    dx = L/(ix-1)
    dy = W/(jx-1)
    K = 2*(1/dx**2 + 1/dy**2) # 1/h^2
    max_iter = 10000

    # Initializing arrays
    x = np.linspace(0.0, L, ix)
    y = np.linspace(0.0, W, jx)
    phi = np.zeros((ix, jx))
    phi_new = np.zeros((ix, jx))
    max_res_list = []

    # Boundary Condition
    phi[:, -1] = np.ones(ix)
    phi_new[:, -1] = np.ones(ix)

    for k in range(max_iter):

        # Gauss-Seidel
        for j in range(1, jx-1):
            for i in range(1, ix-1):
                 phi_new[i, j] = 1/K * ((phi[i+1, j] + phi_new[i-1, j])/dx**2 + (phi[i, j+1] + phi_new[i, j-1])/dy**2)

        # Computing max residual
        max_res = np.max(np.abs(phi_new - phi))
        max_res_list.append(max_res)

        if max_res < res:
            break

        # Shifting computation window back
        phi = phi_new.copy()

    # Plotting
    plt.figure('Gauss-Seidel Max Residual')
    plt.semilogy(np.arange(k+1), max_res_list, label=f"ix = jx = {ix}")
    plt.xlabel('Iterations')
    plt.ylabel('Max Residual')
    plt.legend()
    
    plt.figure('Gauss-Seidel Line Plot')
    idy = np.argmin(np.abs(y - W/2))
    plt.plot(x, phi[:, idy], label=f"ix = jx = {ix}")
    plt.xlabel('x')
    plt.ylabel('phi(x, W/2)')
    if ix == 50:
        plt.plot(x, phi_exact, label="Exact Solution")
    plt.legend()


def sor(ix, jx):

    # Constants
    w = 1.8 # relaxation factor
    res = 1e-5
    dx = L/(ix-1)
    dy = W/(jx-1)
    K = 2*(1/dx**2 + 1/dy**2) # 1/h^2
    max_iter = 10000

    # Initializing arrays
    x = np.linspace(0.0, L, ix)
    y = np.linspace(0.0, W, jx)
    phi = np.zeros((ix, jx))
    phi_temp = np.zeros((ix, jx))
    phi_new = np.zeros((ix, jx))
    max_res_list = []

    # Boundary Condition
    phi[:, -1] = np.ones(ix)
    phi_temp[:, -1] = np.ones(ix)
    phi_new[:, -1] = np.ones(ix)

    for k in range(max_iter):

        # Successive Over-Relaxation
        for j in range(1, jx-1):
            for i in range(1, ix-1):
                 phi_temp[i, j] = 1/K * ((phi[i+1, j] + phi_new[i-1, j])/dx**2 + (phi[i, j+1] + phi_new[i, j-1])/dy**2)
                 phi_new[i, j] = phi[i, j] + w * (phi_temp[i, j] - phi[i, j])

        # Computing max residual
        max_res = np.max(np.abs(phi_new - phi))
        max_res_list.append(max_res)

        if max_res < res:
            break

        # Shifting computation window back
        phi = phi_new.copy()

    # Plotting
    plt.figure('SOR Max Residual')
    plt.semilogy(np.arange(k+1), max_res_list, label=f"ix = jx = {ix}")
    plt.xlabel('Iterations')
    plt.ylabel('Max Residual')
    plt.legend()
    
    plt.figure('SOR Line Plot')
    idy = np.argmin(np.abs(y - W/2))
    plt.plot(x, phi[:, idy], label=f"ix = jx = {ix}")
    plt.xlabel('x')
    plt.ylabel('phi(x, W/2)')
    if ix == 50:
        plt.plot(x, phi_exact, label="Exact Solution")
    plt.legend()


def slor(ix, jx):

    # Constants
    w = 1.8 # relaxation factor
    res = 1e-5
    dx = L/(ix-1)
    dy = W/(jx-1)
    K = 2*(1/dx**2 + 1/dy**2) # 1/h^2
    max_iter = 10000

    # Initializing arrays
    x = np.linspace(0.0, L, ix)
    y = np.linspace(0.0, W, jx)
    phi = np.zeros((ix, jx))
    phi_temp = np.zeros((ix, jx))
    phi_new = np.zeros((ix, jx))
    p = np.ones(jx) / dy**2
    q = np.ones(jx) * -K
    r = np.ones(jx) / dy**2
    s = np.zeros(jx)
    max_res_list = []

    # Boundary Condition
    phi[:, -1] = np.ones(ix)
    phi_temp[:, -1] = np.ones(ix)
    phi_new[:, -1] = np.ones(ix)

    for k in range(max_iter):

        # Sweep in x and solve for each line of constant x using the thomas algorithm (to propagate BCs quickly)
        for i in range(1, ix-1):
            s = - (phi[i+1, :] + phi_new[i-1, :]) / dx**2
            phi_temp[i, :] = thomas_algorithm(p, q, r, s)
            phi_new[i, 1:-1] = phi[i, 1:-1] + w * (phi_temp[i, 1:-1] - phi[i, 1:-1])

        # Computing max residual
        max_res = np.max(np.abs(phi_new - phi))
        max_res_list.append(max_res)

        if max_res < res:
            break

        # Shifting computation window back
        phi = phi_new.copy()

    # Plotting
    plt.figure('SLOR Max Residual')
    plt.semilogy(np.arange(k+1), max_res_list, label=f"ix = jx = {ix}")
    plt.xlabel('Iterations')
    plt.ylabel('Max Residual')
    plt.legend()
    
    plt.figure('SLOR Line Plot')
    idy = np.argmin(np.abs(y - W/2))
    plt.plot(x, phi[:, idy], label=f"ix = jx = {ix}")
    plt.xlabel('x')
    plt.ylabel('phi(x, W/2)')
    if ix == 50:
        plt.plot(x, phi_exact, label="Exact Solution")
    plt.legend()


if __name__ == "__main__":
    
    L = 1.0
    W = 1.0

    # Exact Solution
    phi_exact = exact_sol(50, 50)
    '''
    # Jacobi
    jacobi(50, 50)
    jacobi(100, 100)

    # Gauss-Seidel
    gauss_seidel(50, 50)
    gauss_seidel(100, 100)

    # SOR
    sor(50, 50)
    sor(100, 100)
    '''
    # SLOR
    slor(50, 50)
    slor(100, 100)

    plt.show()