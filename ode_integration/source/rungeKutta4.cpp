#include <rungeKutta4.h>

std::vector<std::pair<double, double>> rungeKutta4Solver(
    std::function<double(double, double)> f, // Input function f' = f(x, y)
    std::pair<double, double> xInit, // Initial conditions
    double xEnd, // Solve until x = xEnd
    double h // Step size
){

    std::vector<std::pair<double, double>> solution; // Stores solution
    solution.push_back(xInit);

    double x = xInit.first;
    double y = xInit.second;

    while (x <= xEnd){
        double k1 = h * f(x, y);
        double k2 = h * f(x + h / 2, y + k1 / 2);
        double k3 = h * f(x + h / 2, y + k2 / 2);
        double k4 = h * f(x + h, y + k3);
        
        y += (k1 + 2 * k2 + 2 * k3 + k4) / 6; // Update y
        x += h; // Update x

        solution.push_back({x, y});
    }

    return solution;
}