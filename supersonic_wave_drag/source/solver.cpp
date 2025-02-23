#include <solver.h>
#include <iostream>

double sliding_window(double M){

    // Domain
    const double x0 = 0.0;
    const double xEnd = 1.0;
    const double y0 = 0.0;
    const double yEnd = 2.0;

    // Step size
    const int jx = 50;
    const double dy = (yEnd - y0) / jx;
    double dx = dy * std::sqrt(M*M - 1);
    const int ix = (xEnd - x0) / dx;

    // Airfoil constants
    const double e_m = 0.1;

    // Initializing solution matrix
    std::vector<std::vector<double>> phi(3, std::vector<double>(jx+1));

    // Initial Conditions
    for (int i = 0; i < 2; ++i) {
        std::fill(phi[i].begin(), phi[i].end(), 0.0);
    }

    // Solving for next step
    double cd = 0;
    for (int i = 1; i < ix + 1; ++i) {
        for (int j = 1; j < jx; ++j) {
            // Solving next values in the domain
            phi[2][j] = - phi[0][j] + phi[1][j+1] + phi[1][j-1];
        }
        // Airfoil gradient
        double df = 2 * e_m * (1 - 2 * (i - 0.5) * dx);

        // Boundary conditions
        phi[2][jx] = phi[2][jx-1];
        phi[2][0] = phi[2][1] - dy * df;

        // Drag Coefficient
        cd += -4 * (phi[2][0] - phi[1][0]) * df;

        // Updating computational window
        phi[0] = phi[1];
        phi[1] = phi[2];
    }

    return cd;
}


double full_domain(double M){

        // Domain
        const double x0 = 0.0;
        const double xEnd = 1.0;
        const double y0 = 0.0;
        const double yEnd = 2.0;
    
        // Step size
        const int jx = 50;
        const double dy = (yEnd - y0) / jx;
        double dx = dy * std::sqrt(M*M - 1);
        const int ix = (xEnd - x0) / dx;
    
        // Airfoil constants
        const double e_m = 0.1;
    
        // Initializing solution matrix
        std::vector<std::vector<double>> phi(ix + 2, std::vector<double>(jx+1));
    
        // Initial Conditions
        for (int i = 0; i < 2; ++i) {
            std::fill(phi[i].begin(), phi[i].end(), 0.0);
        }
    
        // Solving for next step
        double cd = 0;
        for (int i = 1; i < ix + 1; ++i) {
            for (int j = 1; j < jx; ++j) {
                // Solving next values in the domain
                phi[i+1][j] = - phi[i-1][j] + phi[i][j+1] + phi[i][j-1];
            }
            // Airfoil gradient
            double df = 2 * e_m * (1 - 2 * (i - 0.5) * dx);
    
            // Boundary conditions
            phi[i+1][jx] = phi[i+1][jx-1];
            phi[i+1][0] = phi[i+1][1] - dy * df;
    
            // Drag Coefficient
            cd += -4 * (phi[i+1][0] - phi[i][0]) * df;
        }
    
        return cd;
}