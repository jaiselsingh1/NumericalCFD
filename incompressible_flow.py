import numpy as np
import matplotlib.pyplot as plt


def compute_residual(phi, dx, dy):
    
    res = np.zeros_like(phi)

    for j in range(1, res.shape[1] - 1):
        for i in range(1, res.shape[0] - 1):
            res[i, j] = (phi[i+1, j] - 2*phi[i, j] + phi[i-1, j])/dx**2 + (phi[i, j+1] - 2*phi[i, j] + phi[i, j-1])/dy**2

    max_res = np.max(np.abs(res))

    return max_res


def incompressible_flow(dx, dy):

    # Constants
    w = 1.8 # relaxation factor
    tol = 1e-5
    K = 2*(1/dx**2 + 1/dy**2)
    K_airfoil = 2/dx**2 + 1/dy**2
    max_res = np.inf
    max_iter = 1000
    iter = 0

    # Domain
    L = 5
    W = 5
    x = np.arange(-L, L + dx, dx)
    y = np.arange(0, W + dy, dy)

    # Initializing grid
    phi = np.zeros((len(x), len(y)))
    phi_new = np.zeros((len(x), len(y)))

    while (max_res > tol) and (iter < max_iter):

        iter += 1

        for i in range(1, phi.shape[0] - 1):
            for j in range(0, phi.shape[1] - 1):

                # Check if at bottom BC
                if j == 0:
                    # Check if at airfoil
                    if (x[i] > 0) and (x[i] <= 1):
                        f_prime = 0.6 * (0.14845*x[i]**(-0.5) - 0.126 - 0.7032*x[i] + 0.8529*x[i]**2 - 0.406*x[i]**3)
                        phi_new[i, j] = 1/K_airfoil * ((phi[i+1, j] + phi[i-1, j])/dx**2 + (phi[i, j+1]/dy - f_prime)/dy)
                    else:
                        phi_new[i, j] = 1/K * ((phi[i+1, j] + phi[i-1, j])/dx**2 + 2*phi[i, j+1]/dy**2)
                else:
                    # SOR
                    phi_temp = 1/K * ((phi[i+1, j] + phi_new[i-1, j])/dx**2 + (phi[i, j+1] + phi_new[i, j-1])/dy**2)
                    phi_new[i, j] = phi[i, j] + w * (phi_temp - phi[i, j])

        # Computing max residual
        max_res = compute_residual(phi_new, dx, dy)

        # Shifting computation window back
        phi = phi_new.copy()

    if iter == max_iter:
        print(f"Reached max iterations with max residual = {max_res}")
    else:
        print(f"Converged at iteration {iter}")

if __name__ == "__main__":
    incompressible_flow(dx=0.02, dy=0.02)
    print("done")