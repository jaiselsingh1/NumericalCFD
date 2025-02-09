#include <improvedEulerCauchy.h>

std::vector<std::pair<double, double>> improvedEulerCauchySolver(
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

        double y_predictor = y + h * f(x, y);
        double y_corrector = y + h/2 * (f(x, y) + f(x+h, y_predictor));

        y = y_corrector; // Update y
        x += h; // Update x

        solution.push_back({x, y});
    }

    return solution;
}