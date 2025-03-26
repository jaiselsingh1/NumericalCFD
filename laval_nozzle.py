import numpy as np
import matplotlib.pyplot as plt


def laval_nozzle_numerical_solver(u_outlet, ix):

    # Parameters
    eps = 1e-6
    dx = 1/ix
    dt = dx/-u_outlet

    # Laval Nozzle Geometry
    g = np.zeros(ix+1)
    g_prime = np.zeros(ix+1)
    
    for i in range(ix+1):
        x = i*dx
        g[i] = 1 + 4*(x - 0.5)**2
        g_prime[i] = 8*x - 4

    # Initial Conditions
    u = np.zeros(ix+1)
    u_prev = np.zeros(ix+1)

    # Boundary Conditions
    u[0] = u_outlet
    u[-1] = u_outlet

    # Steady-State Solver
    while not (np.linalg.norm(u - u_prev) < eps):

        u_next = np.zeros(ix+1)
        u_next[0] = u_outlet
        u_next[-1] = u_outlet

        for i in range(1, ix):

            u_minus_half = (u[i-1] + u[i])/2
            u_plus_half = (u[i+1] + u[i])/2
            
            RHS = (1-(u[i]**2)/2)*g_prime[i]/g[i]

            if (u_minus_half > 0) and (u_plus_half > 0):
                u_next[i] = u[i] + dt * (RHS - u_minus_half*(u[i] - u[i-1])/dx)
            elif (u_minus_half > 0) and (u_plus_half < 0):
                u_next[i] = u[i] + dt * (RHS - u_minus_half*(u[i] - u[i-1])/dx - u_plus_half*(u[i+1] - u[i])/dx)
            elif (u_minus_half <= 0) and (u_plus_half >= 0):
                u_next[i] = (u[i] + dt*RHS)/(1 + dt/(2*dx)*(u[i+1] - u[i-1]))
            elif (u_minus_half < 0) and (u_plus_half < 0):
                u_next[i] = u[i] + dt * (RHS - u_plus_half*(u[i+1] - u[i])/dx)

        u_prev = u
        u = u_next

    return u


def laval_nozzle_exact_solver(u_outlet, ix):

    # Parameters
    dx = 1/ix

    # Initializing matrices
    u = np.zeros(ix+1)

    # Exact Solution
    for i in range(ix+1):
        x = i * dx
        u[i] = -np.sqrt(2 - (4 - 2*u_outlet**2)/(1 + 4*(x - 0.5)**2))
    
    return u


if __name__ == "__main__":
    ix = 210
    x = np.linspace(0, 1, ix+1)

    u_num_a = laval_nozzle_numerical_solver(-np.sqrt(1.5), ix)
    u_exact_a = laval_nozzle_exact_solver(-np.sqrt(1.5), ix)

    u_num_b = laval_nozzle_numerical_solver(-1., ix)
    u_exact_b = laval_nozzle_exact_solver(-1., ix)
    
    plt.figure()
    plt.plot(x, u_num_a, label='Numerical')
    plt.plot(x, u_exact_a, label="Exact")
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('u_outlet = -sqrt(1.5)')
    plt.legend()

    plt.figure()
    plt.plot(x, u_num_b, label='Numerical')
    plt.plot(x, u_exact_b, label="Exact")
    plt.xlabel('x')
    plt.ylabel('u')
    plt.title('u_outlet = -1.0')
    plt.legend()

    plt.show()