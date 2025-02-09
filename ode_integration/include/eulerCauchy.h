#ifndef EULERCAUCHY_H
#define EULERCAUCHY_H

#include <vector>
#include <functional>

std::vector<std::pair<double, double>> eulerCauchySolver(
    std::function<double(double, double)> f, // Input function
    std::pair<double, double> x0, // Initial conditions
    double xEnd, // Solve until x = xEnd
    double h // Step size
);

#endif