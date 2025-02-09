#include <eulerCauchy.h>

std::vector<std::pair<double, double>> eulerCauchySolver(
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
        
        y += h * f(x, y); // Update y
        x += h; // Update x

        solution.push_back({x, y});
    }

    return solution;
}